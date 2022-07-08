#' Compute community signatures (eigengenes)
#'
#' Computes the community signatures (eigengenes) for an expression matrix
#' given a particular community membership (kME) matrix
#'
#' @param ex matrix of bulk RNA-seq or microarray gene expression data.
#' This should be in log space and greater than or equal to 0.
#' @param membership_matrix, a community membership (kME) matrix with genes as
#' rows and communities as columns.
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
#' @export
compute_eigengene_matrix <- function(ex,
                                     membership_matrix,
                                     cutoff = 5,
                                     pc_flag = TRUE) {
  ex <- ex[apply(as.matrix(ex), 1, stats::sd) != 0, ]
  membership_matrix <- membership_matrix[rownames(membership_matrix) %in% rownames(ex), ]
  m_genes <- rownames(membership_matrix)
  e_genes <- rownames(ex)

  tEigen <- plyr::aaply(as.matrix(membership_matrix), 2, function(k) {
    tInds <- rank(-k) <= cutoff
    tEx <- ex[e_genes %in% m_genes[tInds], ]
    if (all(tEx == 0)) {
      return(rep(0, ncol(tEx)))
    }

    scaled_ave <- apply(t(scale(t(tEx))), 2, mean)
    if (pc_flag == FALSE) {
      return(scaled_ave)
    } else {
      pc1 <- stats::prcomp(t(tEx))$x[, 1]
      if (stats::cor(scaled_ave, pc1) < 0) {
        pc1 <- -pc1
      }
      return(pc1)
    }
  })
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
prolif_names <- c("TPX2","PRC1","BIRC5","CEP55","MELK","KIF4A","CDC20",
                  "MCM10","HJURP","FOXM1","TOP2A","DLGAP5","KIF2C","KIF14",
                  "ASPM","NEK2","CDCA8","CDKN3","NUF2","CDCA3",
                  "CCNA2","CDCA5","CCNB1","ANLN","TTK","KIF20A","CCNB2")


#' Compute Cell Type Enrichments Using panglaoDB Cell Markers
#'
#' compute cell type enrichments using [panglaoDB cell markers](https://academic.oup.com/database/article/doi/10.1093/database/baz046/5427041?login=true)
#' using Fisher test.
#' @param t_memb `community_membership` or `full_community_membership`
#' values from [icwgcna()]
#' @param K cutoff for top community genes to include for computing enrichment.
#' Used in an AND condition with memb_cut.
#' @param memb_cut cutoff as a membership score threshold for determining top
#' community genes for computing enrichment.  Used in an AND condition with K.
#' @param pangDB panglaoDB cell markers database. Default is to read data from t
#' he url [pangDB_link]
#' @param prolif list of proliferation genes to check.
#' Default is to use [prolif_names]
#'
#' @return Returns a list with the following items:
#' * `top_enr` - the most significant cell type from the enrichment scores.
#' * `full_enr` - all panglaoDB cell type enrichment scores for all communities.
#'
#'
#' @export
#'
#' @examples
#'
#'\dontrun{
#' pangDB <- data.table::fread(pangDB_link)
#' compute_panglaoDB_enrichment(tcell_net$community_membership, pangDB = pangDB)}
#'
compute_panglaoDB_enrichment <- function(t_memb,
                                         K = 100,
                                         memb_cut = .65,
                                         pangDB = data.table::fread(pangDB_link),
                                         prolif = prolif_names) {
  if (!any(class(t_memb) %in% c('matrix', 'data.frame'))) {
    stop("t_memb must be a martix or data.frame")
  }
  if (min(t_memb) < -1 || max(t_memb) > 1) {
    stop("t_memb values can't be <-1 or >1")
  }

  colnames(pangDB) <- make.names(colnames(pangDB))
  if (!all(c('cell.type', 'official.gene.symbol') %in% colnames(pangDB))) {
    stop('expecting pangDB variables "cell.type" and "official.gene.symbol" (after making syntactically valid names using make.names() function)')
  }

  if (!any(rownames(t_memb) %in% prolif)) {
    stop('No rownames of "t_memb" are in the provided "prolif" list')
  }

  c.types <- as.vector(stats::na.omit(unique(pangDB$cell.type)))
  enr <- plyr::adply(t_memb,2, function(x) {
    prolif_overlap <- table(
      rank(-x) <= K & x > memb_cut,
      rownames(t_memb) %in% prolif
    )
    ret <- t(plyr::ldply(c.types,function(ct) {
      if (ncol(prolif_overlap) > 1 &
          (prolif_overlap[2,2] > floor(.33 * length(prolif)))) {
        return(1)
      }

      overlap <- table(
        rank(-x) <= K & x > memb_cut,
        rownames(t_memb) %in%
          pangDB$official.gene.symbol[pangDB$cell.type == ct]
      )
      if (ncol(overlap) == 1) {
        return(1)
      }
      ret2 <- stats::fisher.test(overlap,)$p.val
      return(unlist(ret2))
    }))

    colnames(ret) <- c.types
    return(ret)
  })

  rownames(enr) <- enr[,1]
  enr           <- as.data.frame(t(enr[,-1]))

  top_enr <- plyr::ldply(apply(enr,2,function(x){
    if (min(x) > 0.001) {
      return(data.frame(cell_type = NA,p = NA))
    }

    i <- which(x == min(x))

    ret <- data.frame(cell_type = rownames(enr)[i], p = x[i])
    return(ret)
  }))
  names(top_enr)[1] <- "community"

  return(list(top_enr = top_enr,
              full_enr = enr))
}


