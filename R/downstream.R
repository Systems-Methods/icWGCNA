#' Compute community signatures (eigengenes)
#'
#' Computes the community signatures (eigengenes) for an expression matrix
#' given a particular community membership (kME) matrix
#'
#' @param ex matrix of bulk RNA-seq or microarray gene expression data.
#' This should be in log space and greater than or equal to 0.
#' @param membership_matrix a community membership (kME) matrix with genes as
#' rows and communities as columns. Often `community_membership` or
#' `full_community_membership` output from [icwgcna()]
#' @param cutoff number of top genes to use when computing community signatures
#' @param pc_flag indicator. TRUE (default) means to use the 1st principal component
#' (corrected for direction). FALSE uses the mean of scaled and centered top genes.
#'
#' @return A matrix with rows being the community signature and columns being samples
#'
#' @details
#' This can be used to compute community signatures in a new expression dataset.
#' Note, community signatures are not corrected by [icwgcna()]
#' iterations so it will not match signatures output from [icwgcna()]
#' if it is run on the network construction dataset.
#'
#' When using these community signatures for modeling it may be best to include
#' interaction terms or use tree based
#' methods since dependencies are not addressed in this output matrix.
#'
#' @examples
#' \dontrun{
#' library("UCSCXenaTools")
#' luad <- getTCGAdata(
#'   project = "LUAD", mRNASeq = TRUE, mRNASeqType = "normalized",
#'   clinical = FALSE, download = TRUE
#' )
#' ex <- as.matrix(data.table::fread(luad$destfiles), rownames = 1)
#'
#' results <- icwgcna(ex)
#'
#' compute_eigengene_matrix(ex, results$community_membership)
#' }
#'
#' @export
compute_eigengene_matrix <- function(ex,
                                     membership_matrix,
                                     cutoff = 5,
                                     pc_flag = TRUE) {
  if (!all(unlist(lapply(ex, is.numeric)))) {
    stop("all 'ex' columns must be numeric")
  }
  if (min(ex) < 0) {
    stop("all values of ex must be >=0")
  }
  if (max(ex) > 100) {
    warning("some values of ex are >100, strongly indicating ex is not in log space")
  }

  rows_to_use <- rownames(membership_matrix) %in% rownames(ex)
  if (!any(rows_to_use)) {
    stop("No matching rownames in ex and membership_matrix")
  }


  SD <- apply(as.matrix(ex), 1, stats::sd)
  # need to remove any 0 SD genes
  SD_zero_index <- SD == 0
  if (any(SD_zero_index)) {
    message(
      "Removing ", sum(SD_zero_index),
      " genes with a 0 standard deviation"
    )
    ex <- ex[!SD_zero_index, , drop = FALSE]
  }

  membership_matrix <- membership_matrix[rows_to_use, ]
  m_genes <- rownames(membership_matrix)
  e_genes <- rownames(ex)

  tEigen <- plyr::aaply(as.matrix(membership_matrix), 2, function(k) {
    tInds <- rank(-k) <= cutoff
    tEx <- ex[e_genes %in% m_genes[tInds], , drop = FALSE]
    if (all(tEx == 0)) {
      return(rep(0, ncol(tEx)))
    }

    scaled_ave <- apply(t(scale(t(tEx))), 2, mean)
    if (!pc_flag) {
      scaled_ave
    } else {
      pc1 <- stats::prcomp(t(tEx))$x[, 1]
      if (stats::cor(scaled_ave, pc1) < 0) {
        pc1 <- -pc1
      }
      pc1
    }
  })
  tEigen
}


#' Current link used to connect to [panglaoDB cell markers](https://academic.oup.com/database/article/doi/10.1093/database/baz046/5427041?login=true)
#'
#' @details
#' pangDB_link = [https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz](https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz)
#'
#' @examples
#'
#' pangDB_link
#'
#' @export
#'
pangDB_link <- "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz"

#' Current list of proliferation genes to check
#'
#' @examples
#'
#' prolif_names
#'
#' @export
#'
prolif_names <- c(
  "TPX2", "PRC1", "BIRC5", "CEP55", "MELK", "KIF4A", "CDC20",
  "MCM10", "HJURP", "FOXM1", "TOP2A", "DLGAP5", "KIF2C", "KIF14",
  "ASPM", "NEK2", "CDCA8", "CDKN3", "NUF2", "CDCA3",
  "CCNA2", "CDCA5", "CCNB1", "ANLN", "TTK", "KIF20A", "CCNB2"
)


#' Compute Cell Type Enrichments Using panglaoDB Cell Markers
#'
#' compute cell type enrichments using [panglaoDB cell markers](https://academic.oup.com/database/article/doi/10.1093/database/baz046/5427041?login=true)
#' using Fisher test.
#' @param membership_matrix a community membership (kME) matrix with genes as
#' rows and communities as columns. Often `community_membership` or
#' `full_community_membership` output from [icwgcna()]
#' @param K cutoff for top community genes to include for computing enrichment.
#' Used in an AND condition with memb_cut.
#' @param memb_cut cutoff as a membership score threshold for determining top
#' community genes for computing enrichment.  Used in an AND condition with K.
#' @param pangDB panglaoDB cell markers database. Default is to read data from t
#' he url [pangDB_link]
#' @param prolif list of proliferation genes to check.
#' Default is to use [prolif_names]
#' @param p_cut p value cutoff for determining importance. If all p values are
#' below `p_cut` for a community, no cell type is selected
#'
#' @return Returns a list with the following items:
#' * `top_enr` - the most significant cell type from the enrichment scores.
#' * `full_enr` - all panglaoDB cell type enrichment scores for all communities.
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library("UCSCXenaTools")
#' luad <- getTCGAdata(
#'   project = "LUAD", mRNASeq = TRUE, mRNASeqType = "normalized",
#'   clinical = FALSE, download = TRUE
#' )
#' ex <- as.matrix(data.table::fread(luad$destfiles), rownames = 1)
#' results <- icwgcna(ex)
#'
#' pangDB <- data.table::fread(pangDB_link)
#' compute_panglaoDB_enrichment(results$community_membership, pangDB = pangDB)
#' }
#'
compute_panglaoDB_enrichment <- function(membership_matrix,
                                         K = 100,
                                         memb_cut = .65,
                                         pangDB = data.table::fread(pangDB_link),
                                         prolif = prolif_names,
                                         p_cut = 0.001) {
  if (!any(class(membership_matrix) %in% c("matrix", "data.frame"))) {
    stop("membership_matrix must be a martix or data.frame")
  }
  if (min(membership_matrix) < -1 || max(membership_matrix) > 1) {
    stop("membership_matrix values can't be <-1 or >1")
  }

  colnames(pangDB) <- make.names(colnames(pangDB))
  if (!all(c("cell.type", "official.gene.symbol") %in% colnames(pangDB))) {
    stop('expecting pangDB variables "cell.type" and "official.gene.symbol" (after making syntactically valid names using make.names() function)')
  }

  c.types <- as.vector(stats::na.omit(unique(pangDB$cell.type)))
  enr <- plyr::adply(membership_matrix, 2, function(x) {
    prolif_overlap <- table(
      rank(-x) <= K & x > memb_cut,
      rownames(membership_matrix) %in% prolif
    )
    if (nrow(prolif_overlap) < 2 ||
      ncol(prolif_overlap) < 2 ||
      prolif_overlap[2, 2] > floor(.33 * length(prolif))) {
      ret <- t(rep(1, length(c.types)))
    } else {
      ret <- t(plyr::ldply(c.types, function(ct) {
        overlap <- table(
          rank(-x) <= K & x > memb_cut,
          rownames(membership_matrix) %in%
            pangDB$official.gene.symbol[pangDB$cell.type == ct]
        )
        if (ncol(overlap) == 1 || nrow(overlap) == 1) {
          return(1)
        }
        ret2 <- stats::fisher.test(overlap)$p.val
        return(unlist(ret2))
      }))
    }
    colnames(ret) <- c.types
    ret
  })

  rownames(enr) <- enr[, 1]
  enr <- as.data.frame(t(enr[, -1]))

  top_enr <- plyr::ldply(apply(enr, 2, function(x) {
    if (min(x) > p_cut) {
      return(data.frame(cell_type = NA, p = NA))
    }

    i <- which(x == min(x))

    ret <- data.frame(cell_type = rownames(enr)[i], p = x[i])
    return(ret)
  }))
  names(top_enr)[1] <- "community"

  return(list(
    top_enr = top_enr,
    full_enr = enr
  ))
}



#' Compute MSigDB Collection Enrichments for each community
#'
#' Compute MSigDB Collection enrichments using msigdbr
#' using Fisher test.
#' @param membership_matrix a community membership (kME) matrix with genes as
#' rows and communities as columns. Often `community_membership` or
#' `full_community_membership` output from [icwgcna()]
#' @param K cutoff for top community genes to include for computing enrichment.
#' Used in an AND condition with `memb_cut`.
#' @param memb_cut cutoff as a membership score threshold for determining top
#' community genes for computing enrichment.  Used in an AND condition with `K`.
#' @param cats MSigDB collections to use. We recommend only using
#' H, C3, C6, C7, C8 and avoiding C1, C2, C4, C5 for speed
#' @param p_cut p value cutoff for determining importance. If all p values are
#' below `p_cut` for a community, no cell type is selected
#'
#' @return Returns a list with the following items:
#' * `top_enr` - a data.frame of the most significant enrichment for each
#' MSigDB collection.
#' * `full_enr` - all MSigDB enrichment scores for all communities for the
#' selected collections.
#'
#' @details
#' The function can use parallel and distributed processing in R, via the
#' [foreach][foreach::foreach()] package.
#'
#' GSEA Enrichment testing assumes differential analysis so we are using simple
#' fisher tests instead (actually hypergeometric for speed)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library("UCSCXenaTools")
#' luad <- getTCGAdata(
#'   project = "LUAD", mRNASeq = TRUE, mRNASeqType = "normalized",
#'   clinical = FALSE, download = TRUE
#' )
#' ex <- as.matrix(data.table::fread(luad$destfiles), rownames = 1)
#' results <- icwgcna(ex)
#'
#' # Running with whatever parallel processing is already set up
#' compute_MSigDB_enrichment(results$community_membership)
#'
#' # Using doParallel package to set up parallel processing
#' cl <- parallel::makePSOCKcluster(2)
#' doParallel::registerDoParallel(cl)
#' compute_MSigDB_enrichment(results$community_membership)
#'
#' # Using doMC package to set up parallel processing (not available in Windows)
#' doMC::registerDoMC()
#' compute_MSigDB_enrichment(results$community_membership)
#' }
#'
compute_MSigDB_enrichment <- function(membership_matrix,
                                      K = 100,
                                      memb_cut = .65,
                                      cats = c("H", "C3", "C6", "C7", "C8"),
                                      p_cut = 0.001) {
  needed_packages <- c('msigdbr', 'foreach', 'tidyr')
  missing_packages <- !vapply(needed_packages,
                          FUN = requireNamespace, quietly = TRUE,
                          FUN.VALUE = logical(1))
  if (any(missing_packages)) {
    stop('Must have the following R packages installed for this function: ',
         paste0(names(missing_packages[missing_packages]), collapse = ', '))
  }

  if (!any(class(membership_matrix) %in% c("matrix", "data.frame"))) {
    stop("membership_matrix must be a martix or data.frame")
  }
  if (min(membership_matrix) < -1 || max(membership_matrix) > 1) {
    stop("membership_matrix values can't be <-1 or >1")
  }

  m_df_simp <- msigdbr::msigdbr(species = "Homo sapiens")
  gs_cat_levels <- sort(unique(m_df_simp$gs_cat))

  if (all(!cats %in% gs_cat_levels)) {
    stop('No "cats" found in MSigDB. Must use at least one of: ',
         paste0(gs_cat_levels, collapse = ', '))
  }
  if (any(!cats %in% gs_cat_levels)) {
    warning('The following "cats" are not in MSigDB: ',
            paste0(cats[!cats %in% gs_cat_levels], collapse = ', '))
    cats <- cats[cats %in% gs_cat_levels]
  }

  sig_list <- lapply(cats, function(xx) {
    tmp_data <- m_df_simp[m_df_simp$gs_cat == xx, ]
    split(x = tmp_data$gene_symbol, f = tmp_data$gs_name)
  })
  names(sig_list) <- cats

  if (is.null(foreach::getDoParName())) {
    message('No parallel processing has been detected')
     `%d%` <- foreach::`%do%`
  } else {
    message('Using ', foreach::getDoParName(), ' with ',
            foreach::getDoParWorkers(), ' workers')
    `%d%` <- foreach::`%dopar%`
  }

  top_genes <- apply(membership_matrix, 2,
                     function(meta) {
                       rank(-meta) <= K & meta > memb_cut
                     })
  all_genes <- rownames(top_genes)

  overlap_enrichment <- plyr::ldply(
    names(sig_list),
    .fun = function(xx) {
      sig <- sig_list[[xx]]
      message("working on ", xx)

      goGenes <- as.data.frame(lapply(
        sig,
        function(go) {
          all_genes %in% go
        }))

      top_genes <- top_genes

      i <- NA
      enr <- foreach::foreach(
        i = 1:ncol(top_genes),
        .combine = 'rbind') %d% {
          q_h     <- t(as.matrix(goGenes)) %*% top_genes[, i]
          m_h     <- colSums(goGenes)
          n_h     <- nrow(goGenes) - m_h
          k_h     <- rep(sum(top_genes[, i]), length(q_h))
          hyper_p <- stats::phyper(q_h - 1, m_h, n_h, k_h, lower.tail = FALSE)

          data.frame(
            community = colnames(top_genes)[i],
            go = colnames(goGenes),
            top_comm_gene_n = k_h,
            go_n = m_h,
            overlap = q_h,
            p = hyper_p
          )
        }
      enr$cat <- xx
      enr
    })

  cat_cnts <- table(overlap_enrichment$community,
                    overlap_enrichment$cat)
  correct  <- all(apply(cat_cnts,2,function(x) {length(unique(x)) == 1}))
  if (!correct) {
    stop("problem computing enrichments, try a different parallel computing setup")
  }

  best_of_cat <- plyr::ddply(
    overlap_enrichment,
    .variables = c("community","cat"),
    function(x) {
      ind <- which(x$p == min(x$p,na.rm = TRUE));
      if (length(ind) > 1 || min(x$p, na.rm = TRUE) > p_cut) {
        best            <- "";
        p_val           <- NA
        top_comm_gene_n <- NA
        go_n            <- NA
        overlap         <- NA
      } else {
        best <- x$go[ind];
        p_val <- x$p[ind]
        top_comm_gene_n <- x$top_comm_gene_n[ind]
        go_n            <- x$go_n[ind]
        overlap         <- x$overlap[ind]
      }
      data.frame(
        best = best,
        top_comm_gene_n,
        go_n,
        overlap,
        p = p_val
      )
    })

  best <- tidyr::pivot_wider(best_of_cat,
                             id_cols = 'community',
                             names_from = 'cat',
                             values_from = c('best', 'top_comm_gene_n',
                                             'go_n', 'overlap', 'p'))
  # reorder col
  best <- best[, c(1, unlist(lapply(cats,
                             function(xx) {
                               grep(paste0('_', xx), colnames(best))
                             }))
                   )]

  list(top_enr = best, full_enr = overlap_enrichment)
}


#' Compute Cell Type Enrichments Using xCell Cell Markers
#'
#' Compute cell type enrichments using
#' [xCell cell markers](https://github.com/dviraran/xCell) with Fisher test.
#' @param membership_matrix a community membership (kME) matrix with genes as
#' rows and communities as columns. Often `community_membership` or
#' `full_community_membership` output from [icwgcna()]
#' @param K cutoff for top community genes to include for computing enrichment.
#' Used in an AND condition with `memb_cut`.
#' @param memb_cut cutoff as a membership score threshold for determining top
#' community genes for computing enrichment. Used in an AND condition with `K`.
#' @param p_cut p value cutoff for determining importance. If all p values are
#' below `p_cut` for a community, no cell type is selected
#'
#' @return Returns a list with the following items:
#' * `top_enr` - the most significant cell type from the enrichment scores.
#' * `full_enr` - all panglaoDB cell type enrichment scores for all communities.
#'
#' @details note that this is distinct from running xCell which is run on expression data. This
#' is an enrichment using their cell markers.
#' @export
#'
#' @examples
#' \dontrun{
#' library("UCSCXenaTools")
#' luad <- getTCGAdata(
#'   project = "LUAD", mRNASeq = TRUE, mRNASeqType = "normalized",
#'   clinical = FALSE, download = TRUE
#' )
#' ex <- as.matrix(data.table::fread(luad$destfiles), rownames = 1)
#' results <- icwgcna(ex)
#'
#' compute_xCell_enrichment(results$community_membership)
#' }
#'
compute_xCell_enrichment <- function(membership_matrix,
                                     K = 100,
                                     memb_cut = .65,
                                     p_cut = 0.001) {
  needed_packages <- c('xCell')
  missing_packages <- !vapply(needed_packages,
                              FUN = requireNamespace, quietly = TRUE,
                              FUN.VALUE = logical(1))
  if (any(missing_packages)) {
    stop('Must have the following R packages installed for this function: ',
         paste0(names(missing_packages[missing_packages]), collapse = ', '))
  }


  if (!any(class(membership_matrix) %in% c("matrix", "data.frame"))) {
    stop("membership_matrix must be a martix or data.frame")
  }
  if (min(membership_matrix) < -1 || max(membership_matrix) > 1) {
    stop("membership_matrix values can't be <-1 or >1")
  }

  markers <- xCell::xCell.data$signatures
  tNames  <- names(markers)
  markers <- plyr::llply(names(markers),function(n){markers[[n]]@geneIds})
  names(markers) <- tNames

  c_names_long <- gsub("%.*$","",names(markers))
  c_names <- unique(c_names_long)
  c_types <- plyr::llply(c_names,
                         function(x) {
                           tmp_index <- which(x == c_names_long)
                           unique(as.vector(unlist(markers[tmp_index])))
                         })
  names(c_types) <- c_names

  enr <- plyr::adply(membership_matrix, 2, function(x) {
    ret <- t(plyr::ldply(names(c_types), function(ct) {
      overlap <- table(
        rank(-x) <= K & x > memb_cut,
        rownames(membership_matrix) %in%
          c_types[[ct]]
      )
      if (ncol(overlap) == 1 || nrow(overlap) == 1) {
        return(1)
      }
      ret2 <- stats::fisher.test(overlap)$p.val
      return(unlist(ret2))
    }))
    colnames(ret) <- names(c_types)
    ret
  })

  rownames(enr) <- enr[, 1]
  enr <- as.data.frame(t(enr[, -1]))

  top_enr <- plyr::ldply(apply(enr, 2, function(x) {
    if (min(x) > p_cut) {
      return(data.frame(cell_type = NA, p = NA))
    }

    i <- which(x == min(x))

    ret <- data.frame(cell_type = rownames(enr)[i], p = x[i])
    return(ret)
  }))
  names(top_enr)[1] <- "community"

  return(list(
    top_enr = top_enr,
    full_enr = enr
  ))
}


#' Display UMAP of Community Membership with text overlays
#'
#' @param membership_matrix a community membership (kME) matrix with genes as
#' rows and communities as columns. Often `community_membership` or
#' `full_community_membership` output from [icwgcna()]
#' @param community_memb_cut_main main cutoff of a membership score threshold
#' for filtering out communities to include in the plot. Any communities with
#' less than `community_n_main` number of genes greater than
#' `community_memb_cut_main` are filtered out.
#' @param community_n_main the number of genes that must have membership scores
#'  great than `community_memb_cut_main` in order for a community to be
#' kept in the plot.
#' @param community_memb_cut_secondary cutoff of a membership score threshold
#' for filtering genes to include in the plot. Any communities with
#' less than `community_n_secondary` number of genes greater than
#' `community_memb_cut_secondary` are filtered out.
#' @param community_n_secondary the number of genes that must have membership
#' scores great than `community_memb_cut_secondary` in order for a community
#' to be kept in the plot.
#' @param gene_memb_cut_main main cutoff of a membership score threshold for
#' filtering genes to include in the plot. Any genes with an absolute membership
#' score greater than `gene_memb_cut_main` in any community is included
#' @param gene_memb_cut_secondary secondary cutoff of a membership score
#' threshold for filtering genes to include in the plot. Any genes with an
#' absolute membership score greater than `gene_memb_cut_main` in more than 1
#' community is included
#' @param community_labels a data.frame with the text to display over gene
#' communities. Expects first column to match column names of `membership_matrix`
#' and second column to contain text of labels associated with each community
#' @param umap_specs configuration for UMAP (default is [umap::umap.defaults]).
#' To use custom specs copy [umap::umap.defaults] and make specific changes
#' (see `Examples`)
#'
#' @return Returns a list with the following items:
#' * `umap_w_annotation` - a UMAP plot of genes with labeled clusters overlaid
#' * `umap_w_legend` - a UMAP plot of genes with a legend and no overlaid labels
#' * `layout` - the UMAP layout to enable customized user plotting
#'
#' @details
#'
#' # Filtering
#'
#' In order to facilitate viewing, communities and genes are filtered prior to
#' plotting. Gene filtering is performed after community filtering.
#'
#' For inclusion either the main or secondary parameters can be satisfied. For
#' example, using the defaults here a community must have either 20 genes over
#' 0.7 kME or 5 genes over 0.8 kME to be included. After community filtering,
#' gene filtering is performed. Again using the defaults here, a gene must
#' have at least one community over 0.75 kME or two communities over 0.65 kME
#' to be included.
#'
#' Setting `community_memb_cut_main` and `gene_memb_cut_main` to 0 will force
#' all communities and genes to be included in the UMAP
#'
#' # Output
#'
#' Both `umap_w_annotation` and `umap_w_legend` outputs are
#' [ggplot][ggplot2::ggplot()] objects, so plot details can easily be added or
#' modified (i.e. `umap_results$umap_w_legend + theme_classic()`).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library("UCSCXenaTools")
#' luad <- getTCGAdata(
#'   project = "LUAD", mRNASeq = TRUE, mRNASeqType = "normalized",
#'   clinical = FALSE, download = TRUE
#' )
#' ex <- as.matrix(data.table::fread(luad$destfiles), rownames = 1)
#' results <- icwgcna(ex)
#'
#' umap_results <- make_network_umap(results$community_membership)
#' umap_results$umap_w_annotation
#' umap_results$umap_w_legend + theme(legend.position = 'top')
#'
#' # can adjust umap specifications
#' custom_umap_specs <- umap::umap.defaults
#' custom_umap_specs$n_neighbors <- 20
#'
#' umap_results <- make_network_umap(results$community_membership,
#'                                   umap_specs = custom_umap_specs)
#' }
#'
make_network_umap <- function(membership_matrix,
                              community_memb_cut_main = 0.7,
                              community_n_main = 20,
                              community_memb_cut_secondary = 0.8,
                              community_n_secondary = 5,
                              gene_memb_cut_main = 0.75,
                              gene_memb_cut_secondary = 0.65,
                              community_labels = NULL,
                              umap_specs = umap::umap.defaults) {

  needed_packages <- c('ggplot2', 'rlang', 'umap')
  missing_packages <- !vapply(needed_packages,
                              FUN = requireNamespace, quietly = TRUE,
                              FUN.VALUE = logical(1))
  if (any(missing_packages)) {
    stop('Must have the following R packages installed for this function: ',
         paste0(names(missing_packages[missing_packages]), collapse = ', '))
  }

  if (!any(class(membership_matrix) %in% c("matrix", "data.frame"))) {
    stop("membership_matrix must be a martix or data.frame")
  }
  if (min(membership_matrix) < -1 || max(membership_matrix) > 1) {
    stop("membership_matrix values can't be <-1 or >1")
  }

  if (!is.null(community_labels)) {
    if (!any(class(community_labels) == 'data.frame') ||
        ncol(community_labels) != 2 ||
        !any(colnames(community_labels) == 'community')) {
      stop('community_labels must be a data.frame with 2 columns and a column named "community"')
    }
  }

  col_inds <- (apply(abs(membership_matrix) > community_memb_cut_main, 2, sum) >=
                   community_n_main) |
    (apply(abs(membership_matrix) > community_memb_cut_secondary, 2, sum) >=
       community_n_secondary)
  if (sum(col_inds) < 2) {
    stop('Must have at least 2 communities after filtering. Try less restrictive cutoffs.')
  }
  kME_mat    <- membership_matrix[, col_inds]

  row_inds   <- apply(kME_mat > gene_memb_cut_main, 1, any) |
    (apply(kME_mat > gene_memb_cut_secondary, 1, sum) > 1)
  if (sum(row_inds) < 2) {
    stop('Must have at least 2 genes after filtering. Try less restrictive cutoffs.')
  }
  kME_mat    <- kME_mat[row_inds, ]

  message("Filtering from ", ncol(membership_matrix),
          " communites to ", ncol(kME_mat), " communities for plotting.")
  message("Then filtering from ", nrow(membership_matrix),
          " genes to ", nrow(kME_mat), " genes for plotting.")

  memb_u     <- umap::umap(kME_mat, config = umap_specs)
  layout_df  <- as.data.frame(memb_u$layout)
  names(layout_df) <- c("UMAP1","UMAP2")
  Community  <- colnames(kME_mat)[apply(kME_mat,1,
                                        function(x) {which(x == max(x))})]
  if (!is.null(community_labels)) {
    add_lab <- community_labels[match(Community, community_labels$community), 2]
    Community[!is.na(add_lab)]  <- paste(Community[!is.na(add_lab)],
                                         add_lab[!is.na(add_lab)])
  }

  layout_df  <- cbind(layout_df, Community = Community)

  u_plot <- ggplot2::ggplot(layout_df,
                            ggplot2::aes(x = !!rlang::sym('UMAP1'),
                                         y = !!rlang::sym('UMAP2'),
                                         color = !!rlang::sym('Community'))) +
    ggplot2::geom_point(size = .75) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

  cell_type_locs <- plyr::ddply(layout_df, "Community",
                                function(x) {
                                  data.frame(UMAP1 = stats::median(x$UMAP1),
                                             UMAP2 = stats::median(x$UMAP2))
                                })
  cell_type_locs <- cell_type_locs[order(cell_type_locs$UMAP1,
                                         cell_type_locs$UMAP2,
                                         decreasing = T),]
  labeled_u_plot <- ggplot2::ggplot(
    layout_df,
    ggplot2::aes(x = !!rlang::sym('UMAP1'),
                 y = !!rlang::sym('UMAP2'),
                 color = !!rlang::sym('Community'))
    ) +
    ggplot2::geom_point(size = .75) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 5),
                   legend.position = "none") +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(size = 5)
    )) +
    ggplot2::annotate("label",x = cell_type_locs$UMAP1,
                      y = cell_type_locs$UMAP2,
                      label = cell_type_locs$Community,
                      size = 2,
                      fill = "white")

  list(umap_w_annotation = labeled_u_plot,
       umap_w_legend = u_plot,
       layout = layout_df)
}



#' Identify Top Gene of Communities that are unique (only belong to one community)
#'
#' @param membership_matrix a community membership (kME) matrix with genes as
#' rows and communities as columns. Often `community_membership` or
#' `full_community_membership` output from [icwgcna()]
#' @param K number of unique genes to find.
#' @param maxIt maximum number of iterations to avoid looping to meaningless top genes
#'
#' @return data.frame with key = community and value = top gene ID. Top genes will only
#' be associated with one module.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' unique_top_genes <- find_unique_top_genes(
#'     results$community_membership,
#'     K = 100,
#'     maxIt = 10)
#' }
#'
find_unique_top_genes <- function(membership_matrix,
                                  K = 10,
                                  maxIt = 10) {
  if (!any(class(membership_matrix) %in% c("matrix", "data.frame"))) {
    stop("membership_matrix must be a martix or data.frame")
  }
  if (min(membership_matrix) < -1 || max(membership_matrix) > 1) {
    stop("membership_matrix values can't be <-1 or >1")
  }

  converged <- FALSE
  it        <- 1
  while (it <= maxIt & !converged) {
    metaGeneRanks   <- t(plyr::aaply(-as.matrix(membership_matrix),2,rank))
    multiMembInd    <- apply(metaGeneRanks <= K,1,sum) >  1
    nFound <- sum(apply(metaGeneRanks <= K, 1,any) & !multiMembInd)
    message(nFound," unique genes found with ",
            sum(multiMembInd)," non-unique")
    if (sum(multiMembInd) == 0) {converged <- TRUE}
    membership_matrix <- membership_matrix[!multiMembInd,]
    it <- it + 1
  }
  if (nrow(membership_matrix) == 0) {
    stop('No unique genes found. Try reducing "K"')
  }
  metaGeneRanks     <- t(plyr::aaply(-as.matrix(membership_matrix),2,rank))
  members           <- plyr::ldply(plyr::alply(metaGeneRanks <= K, 2, which),
                                   names)
  rownames(members) <- members[,1]
  members           <- as.data.frame(t(members[,-1]))
  temp              <- tidyr::gather(members)
  ret               <- plyr::ddply(temp, .variables = "value",
                                   .fun = function(x){paste(sort(x$key),
                                                            collapse = ";")})
  ret <- data.frame(key = ret[,2], value = ret[,1])
  ret <- ret[order(ret$key),]
  ret
}


#' Map Eigengenes on Seurat Object
#'
#' Use scores to calculate module scores for feature expression
#' programs in single cells and applies to Seurat object using
#' [UCell::AddModuleScore_UCell()]
#'
#' @param seurat_obj Seurat Object
#' @param membership a data.frame or matrix with genes as
#' rows and communities as columns. Often `community_membership` or
#' `full_community_membership` output from [icwgcna()]. Doesn't have to be
#' kME scores.
#' @param cutoff_method should cutoff be based on a value, number of
#' top genes, or either method, whichever yields a smaller gene subset
#' @param value_cutoff value cutoff (ignored when `cutoff_method` = "top_gene")
#' @param top_genes_cutoff number of top genes to include
#' (ignored when `cutoff_method` = "value")
#' @param assay Seurat object assay element to use
#' @param slot Pull out data from this slot of the Seurat object
#' @param prefix prefix to add to column names of the Seurat object meta.data
#' @param suffix suffix to add to column names of the Seurat object meta.data
#' @param ncores Number of processors to parallelize computation for
#' [UCell::AddModuleScore_UCell()]
#'
#' @return
#' Seurat Object with additional meta.data columns of community results
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' unique_top_genes <- map_eigengenes_on_seurat(
#'     seurat_obj
#'     results$community_membership
#'     )
#' }
map_eigengenes_on_seurat <- function(
    seurat_obj,
    membership,
    cutoff_method = c('value', 'top_gene', 'both'),
    value_cutoff = 0.75,
    top_genes_cutoff = 10,
    assay = "RNA",
    slot = "data",
    prefix = NULL,
    suffix = "_UCell",
    ncores = 1){

  cutoff_method <- match.arg(cutoff_method)

  needed_packages <- c("Seurat", "UCell")
  missing_packages <- !vapply(needed_packages,
                              FUN = requireNamespace, quietly = TRUE,
                              FUN.VALUE = logical(1))
  if (any(missing_packages)) {
    stop('Must have the following R packages installed for this function: ',
         paste0(names(missing_packages[missing_packages]), collapse = ', '))
  }

  membership <- as.data.frame(membership[
    row.names(membership) %in% row.names(seurat_obj@assays[[assay]]),
  ])

  message('Collecting signatures...')

  if (cutoff_method %in% c('kme', 'both')) {

    genelist <- lapply(
      membership,
      \(x) {
        reordered_genes <- order(x, decreasing = TRUE)
        reordered_genes_values <- x[reordered_genes]
        reordered_genes <- reordered_genes[
          which(reordered_genes_values >= value_cutoff)
        ]
        res <- row.names(membership)[reordered_genes]
        if (cutoff_method == 'either') {
          if (length(res) > top_genes_cutoff) {
            res <- res[1:top_genes_cutoff]
          }
        }
        return(res)
      })
  } else {
    genelist <- lapply(
      membership,
      \(x) row.names(membership)[order(x, decreasing = TRUE)[1:top_genes_cutoff]])
  }


  genelist <- genelist[!sapply(genelist, function(x) length(x) <= 1)]
  comm_names <- names(genelist)
  if (!is.null(prefix)) {
    comm_names <- paste0(prefix, '_', comm_names)
  }
  names(genelist) <- comm_names

  message('Mapping signatures...')
  UCell::AddModuleScore_UCell(
    obj = seurat_obj,
    features = genelist,
    assay = assay,
    slot = slot,
    name = suffix,
    ncores = ncores
  )
}


