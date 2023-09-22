#' Gene ID Mapping of Expression Dataset
#'
#' Converts between Gene IDs For a given expression dataset and mapping file.
#'
#' @param exprs_data continuous value data.frame or matrix with Gene IDs as
#' rownames
#' @param mapping_file mapping file data.frame with the first column the Gene
#' IDs in rownames of `exprs_data` and the second column the IDs to map to.
#' Can have multiple to multiple linking.
#' @param compress_fun the compression method to use in cases where multiple
#' gene IDs link to a single new gene ID. "highest_mean" and "highest_median"
#' pick the highest row based on mean or median, respectively, while "mean",
#' "median", and "sum" aggregate the duplicate rows based on method chosen.
#' "pc1" and "pc1_unscaled" are using the first principal component in
#' [stats::prcomp()], with `scale.` set to TRUE or FALSE, respectively.
#' @param compress_trans the transformation used when compressing the
#' duplicate rows. For example, "log_exp" would take the log of the data,
#' apply the
#' `compress_fun` method, then transpose back using exp (i.e. geometric mean).
#' @param verbose should messages about the compression process be displayed
#'
#' @details
#' Gene IDs in the `exprs_data` that do not link to the first column of
#' `mapping_file` will be excluded from the final output.
#'
#' When pc1 or pc1_unscaled `compress_fun` are specified the 1st principal
#' component is flipped if there is negative correlation with the duplicate
#' rows.
#'
#'
#' @return
#' a data.frame at the new gene ID level, with compression of duplicate rows
#' as outlined in the `compress_fun` and `compress_trans` parameters
#' @export
#'
#' @examples
#' \dontrun{
#' geo <- GEOquery::getGEO("GSE14333")
#' exprs_data <- geo$GSE14333_series_matrix.txt.gz@assayData$exprs
#'
#' library(hgu133plus2.db)
#' mapping_file     <- as.data.frame(hgu133plus2SYMBOL)
#'
#' gene_mapping(exprs_data, mapping_file, "highest_mean", "none")
#'
#' ### Using a more complicated case for GSE83834
#' temp_dir <- tempdir()
#' geo <- GEOquery::getGEO("GSE83834", destdir = temp_dir)
#' sup_mat <- GEOquery::getGEOSuppFiles("GSE83834", baseDir = temp_dir)
#' exprs_data_raw <- as.data.frame(readxl::read_excel(rownames(sup_mat)))
#' table(table(exprs_data_raw$ID))
#'
#' exprs_data <- exprs_data_raw[,-1]
#' rownames(exprs_data) <- sub('\\.[0-9]*$', '', exprs_data_raw$ID)
#'
#' library('biomaRt')
#' mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#' genes <-  rownames(exprs_data)
#' mapping_file <- getBM(filters= "ensembl_gene_id",
#'                       attributes= c("ensembl_gene_id","hgnc_symbol"),
#'                       values = genes, mart= mart)
#'
#' exprs_data_symbols <- gene_mapping(exprs_data, mapping_file)
#'
#' }
gene_mapping <- function(exprs_data, mapping_file,
                         compress_fun = c("mean", "median", "sum",
                                          "pc1", "pc1_unscaled",
                                          "highest_mean","highest_median"),
                         compress_trans = c("none", "log_exp", "exp_log"),
                         verbose = TRUE
                         ) {
  compress_fun <- match.arg(compress_fun)
  compress_trans <- match.arg(compress_trans)

  mapping_file <- stats::na.omit(mapping_file[, 1:2])
  mapping_file <- mapping_file[mapping_file[,2] != '', ]

  mapping_file_linked <- merge(
    data.frame(id = rownames(exprs_data)),
    mapping_file,
    by.x = "id",
    by.y = colnames(mapping_file)[1]
  )
  if (nrow(mapping_file_linked) == 0) {
    stop("Could not link rownames(exprs_data) to mapping_file[,1]!")
  }
  colnames(mapping_file_linked)[2] <- "symbol"

  n_linked <- length(unique(mapping_file_linked$id))
  n_symbols <- length(unique(mapping_file_linked$symbol))

  if (verbose) {
    message(nrow(exprs_data) - n_linked,
            ' exprs_data rows could not be linked using the mapping file, ',
            "resulting in ", n_linked, " rows. \nThese rows link to ",
            n_symbols,
            " distinct gene symbols. \nWill compress duplicate rows using the ",
            compress_fun, " method with ",
            ifelse(compress_trans == "none", "no", compress_trans),
            " transformation."
            )
  }

  exprs_data <- exprs_data[rownames(exprs_data) %in% mapping_file_linked$id, ]
  if (compress_trans == "log_exp" && any(exprs_data <= 0)) {
    stop("Can't do log transformation with exprs_data values <= 0")
  }

  # compressing dups
  if (any(duplicated(mapping_file_linked$symbol))) {
    dup_symbols <- unique(
      mapping_file_linked$symbol[duplicated(mapping_file_linked$symbol)]
    )

    dup_compressed <- lapply(dup_symbols, function(xx) {
      tmp_mapping <- mapping_file_linked[mapping_file_linked$symbol == xx,]
      tmp_subset <- exprs_data[tmp_mapping$id,]

      # transforming data for compression
      tmp_subset <- switch(
        compress_trans,
        none = tmp_subset,
        log_exp = log(tmp_subset),
        exp_log = exp(tmp_subset),
      )

      if (compress_fun %in% c("pc1", "pc1_unscaled")) {
        tmp_scale. <- switch(compress_fun, pc1 = TRUE, pc1_unscaled = FALSE)
        tmp_ave <- switch(
          compress_fun,
          pc1 = apply(t(scale(t(tmp_subset))), 2, mean),
          pc1_unscaled = apply(tmp_subset, 2, mean)
        )
        tmp_compressed <- stats::prcomp(t(tmp_subset),
                                        scale. = tmp_scale.)$x[, 1]
        if (stats::cor(tmp_ave, tmp_compressed) < 0) {
          tmp_compressed <- -tmp_compressed
        }
      }else if (compress_fun %in% c("highest_mean", "highest_median")) {
        tmp_fun <- substr(compress_fun,
                          regexpr('_', compress_fun) + 1,
                          nchar(compress_fun))
        tmp_stat <- apply(tmp_subset, 1, tmp_fun)
        tmp_compressed <- tmp_subset[
          tmp_stat == max(tmp_stat), , drop = FALSE
        ][1, ]
      } else {
        tmp_compressed <- apply(tmp_subset, 2, compress_fun)
      }
      # need to convert back
      tmp_compressed <- switch(
        compress_trans,
        none = tmp_compressed,
        log_exp = exp(tmp_compressed),
        exp_log = log(tmp_compressed),
      )
    })
    dup_df <- do.call(rbind.data.frame, dup_compressed)
    rownames(dup_df) <- dup_symbols
    colnames(dup_df) <- colnames(exprs_data)

    # getting non-duplicate cases
    non_dup_symbols <- setdiff(mapping_file_linked$symbol, dup_symbols)
    non_dup_ids <- mapping_file_linked$id[match(non_dup_symbols,
                                                mapping_file_linked$symbol)
      ]
    single_df <- exprs_data[match(non_dup_ids, rownames(exprs_data)), ]
    rownames(single_df) <- non_dup_symbols

    # combining dup and non-dup and sorting
    final_data <- rbind(dup_df, single_df)
    final_data <- final_data[order(rownames(final_data)), ]

  } else {
    final_data <- exprs_data[match(mapping_file_linked$id,
                                   rownames(exprs_data)), ]
    rownames(final_data) <- mapping_file_linked$symbol
  }

  final_data
}
