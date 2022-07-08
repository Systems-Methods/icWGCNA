
####################################################################################################
# icWGCNA network construction related functions
####################################################################################################

# simple function computing angular distance to use as a adjacency measure.
angularDist <- function(x) {
  # suppress warnings for NaN values on the diagonal (set to 0 later regardless)
  aDist <- suppressWarnings(asin(x)) / (pi / 2)
  return(aDist)
}

# Wrapper of Rfast cora function to include spearman method which is used in calculating adjacency.
# x is a gene expression matrix with rows = genes and columns = samples.
RfastCor_wrapper <- function(x,
                             Method = c("pearson", "spearman")) {
  Method <- match.arg(Method)

  gene_names <- rownames(x)
  x <- as.matrix(x)
  if (Method == "spearman") {
    x <- Rfast::rowRanks(x)
  }

  Cor <- Rfast::cora(Rfast::transpose(x), large = TRUE)

  rownames(Cor) <- colnames(Cor) <- gene_names

  return(Cor)
}


#' RfastTOMdist
#'
#' @param A A N x N adjacency matrix where N is the number of genes. values range from -1:1 with higher values indicating highly similar genes. Often correlation ^exponent, but could be angular distance, mutual information or Euclidean distance
#' @param mat_mult_method method for large matrix multiplication, "Rfast" (default) or "RcppEigen" (see `details` in [`icwgcna()`])
#'
#' @return A N x N distance matrix with smaller values indicating more related genes.
#' @export
#'
#' @description distance based on the topological overlap map from [Ravasz, E., Somera, A., Mongru, D., Oltvai, Z. and Barab´asi, A. (2002). Science](https://pubmed.ncbi.nlm.nih.gov/12202830/)
#' Implemented to using the [Rfast][Rfast::Rfast-package] functions to  speed things up since we will be computing this up to 25 times
#'
#' @examples
RfastTOMdist <- function(A,
                         mat_mult_method = c('Rfast', 'RcppEigen')) {
  mat_mult_method <- match.arg(mat_mult_method)

  diag(A) <- 0
  A[is.na(A)] <- 0
  kk <- Rfast::colsums(A)
  denomHelp <- matrix(kk, ncol = length(kk), nrow = length(kk))
  denomTOM <- Rfast::Pmin(denomHelp, Rfast::transpose(denomHelp)) + (1 - A)
  rm(denomHelp, kk)

  if (mat_mult_method == 'Rfast') {
    numTOM <- Rfast::mat.mult(A, A) + A
  } else {
    numTOM <- eigenMapMatMult(A, A) + A
  }


  out <- 1 - as.matrix(numTOM / denomTOM)
  diag(out) <- 0
  return(out)
}


#
#' fastTOMwrapper
#'
#' @param X a gene expression matrix w each column being one sample and N rows representing genes. X should be in log space (usually between 0 and 20)
#' @param expo the power to raise the similarity measure to default = 6. If set to NULL, angular distance is used to applied to the similarity measure ( asin(x) / (pi/2) ).
#' @param Method "pearson" or "spearman" the similarity measure to use
#' @param mat_mult_method method for large matrix multiplication, "Rfast" (default) or "RcppEigen" (see `details` in [`icwgcna()`])
#'
#' @return A N x N distance matrix with smaller values indicating more related genes.
#' @export
#' @description Topological overlap map (TOM, from [Ravasz, et al (2002). Science](https://pubmed.ncbi.nlm.nih.gov/12202830/)) wrapper using the [Rfast][Rfast::Rfast-package]
#' to speed up calculations. This function can compute TOM for 24K gene matrix in 8 minute of an AWS-EC2 c5.18xlarge instance, though in practice we run it on subsets of most variable genes.
#' Adjacently measure options are Pearson or Spearman correlation raised to a power and angular distance. Note that we only use signed and weighted adjacencies. If users are interested in genes negatively associated
#' with a community they should check community memberships (kME) output from the main function icwgcna().
#'
#'
#' @examples
fastTOMwrapper <- function(X,
                           expo = 6,
                           Method = c("pearson", "spearman"),
                           mat_mult_method = c('Rfast', 'RcppEigen')) {
  Method <- match.arg(Method)
  mat_mult_method <- match.arg(mat_mult_method)

  # compute weighted adjacency matrix using Rfast package
  X <- RfastCor_wrapper(X, Method = Method)
  X[is.na(X)] <- 0
  X[X < 0] <- 0

  if (is.null(expo)) {
    A <- angularDist(X)
  } else {
    A <- X^expo
    A[X < 0] <- 0
  }
  rm(X)

  tom <- RfastTOMdist(A, mat_mult_method = mat_mult_method)
  return(tom)
}


# computes the eigengene for a single community
calcEigenGene <- function(tEx) {
  if (nrow(tEx) == 1) {
    return(unlist(tEx))
  }
  tEx <- as.matrix(tEx)
  pc1 <- stats::prcomp(t(tEx))$x[, 1]
  if (sum(stats::cor(t(tEx), pc1) < 0, na.rm = TRUE) / nrow(tEx) > .5) {
    pc1 <- -pc1
  }
  return(pc1)
}


# wrapper of [dynamicTreeCut::cutreeHybrid] from Langfelder and Horvath's dynamic tree cutting package.
# Instead of merging modules like they do, we'll be dropping correlated modules based on kME kurtosis
# i.e. we'll be selecting between two correlated modules
cutreeHybridWrapper <- function(d,
                                quantCut = 0.75) {
  dend <- stats::hclust(stats::as.dist(d), method = "average")
  refHeight <- stats::quantile(dend$height, .05, type = 1)
  cutHeight <- as.numeric(quantCut * (max(dend$height) - refHeight) + refHeight)
  modules <- dynamicTreeCut::cutreeHybrid(dendro = dend,
                                          distM = d,
                                          cutHeight = cutHeight,
                                          minClusterSize = 5,
                                          pamStage = TRUE,
                                          pamRespectsDendro = FALSE,
                                          verbose = 0)
  return(modules)
}

# drop correlated modules/communities based on kME kurtosis (keeping the community with higher kurtosis)
# Kurts is the metric used for selecting which community to drop
# it can be the size of each community, the average kME of each top X gene in each community or kurtosis each communities kME
# we use kurtosis
dropModuels <- function(eigenGenes,
                        Kurts = NULL,
                        corCut = .95,
                        verbose = FALSE) {
  eigen_Cors <- stats::cor(t(eigenGenes))
  diag(eigen_Cors) <- 0
  while (any(eigen_Cors > corCut)) {
    doubleBreak <- FALSE

    for (i in 1:(nrow(eigen_Cors)))
    {
      for (j in (i + 1):nrow(eigen_Cors))
      {
        if (eigen_Cors[i, j] > corCut) {
          if (verbose) {
            message(paste("threshold", corCut, "exceeded for mods",
                        rownames(eigenGenes)[i], ",",
                        rownames(eigenGenes)[j]))
          }
          if (is.null(Kurts)) {
            removeInd <- j
          } else {
            if (Kurts[i] < Kurts[j]) {
              removeInd <- i
            } else {
              removeInd <- j
            }
            Kurts <- Kurts[-removeInd]
          }

          eigenGenes <- eigenGenes[-removeInd, ]
          eigen_Cors <- stats::cor(t(eigenGenes))
          diag(eigen_Cors) <- 0
          doubleBreak <- TRUE
          break()
        }
      }
      if (doubleBreak) {
        break()
      }
    }
  }
  eigen_Cors <- stats::cor(t(eigenGenes))
  diag(eigen_Cors) <- 0
  message(paste("eigegenes trimmed to", nrow(eigenGenes),
              "due to correlation >", corCut,
              "max eigenCor =", max(signif(eigen_Cors, 2))))
  return(eigenGenes)
}


# This is the simplified, Rfast-based WGCNA function that is used in each iteration
# Rfast is used to speed up both the adjacency calculation and the TOM computation
simpWGCNAsubNet <- function(tEx,
                            expo = 6,
                            Method = c("pearson", "spearman"),
                            n = 15,
                            minMods = 5,
                            corCut = .6,
                            mat_mult_method = c('Rfast', 'RcppEigen')) {
  Method <- match.arg(Method)
  mat_mult_method <- match.arg(mat_mult_method)

  message(paste("Computing", nrow(tEx),
              "x", nrow(tEx),
              "TOM distance for subset of genes with higher variance"))

  TOMd <- fastTOMwrapper(tEx,
                         expo = expo,
                         Method = Method,
                         mat_mult_method = mat_mult_method)
  mods <- cutreeHybridWrapper(TOMd)$labels
  modSz <- table(mods)
  if (sum(modSz >= n) < minMods) {
    retMods <- names(modSz[rank(-modSz) <= minMods])
  } else {
    retMods <- names(modSz[modSz >= n])
  }

  # omit the catch all module for genes that do not belong to a community
  retMods <- retMods[retMods != "0"]
  message(paste("number of modules found is", length(retMods)))

  if (length(retMods) == 0) {
    eigenGenes <- NULL
  } else {
    eigenGenes <- as.matrix(plyr::ldply(.data = retMods,
                                         .fun = function(m) {
                                         eg <- calcEigenGene(tEx[mods == m, ])
                                         return(eg)
                                       }))
    rownames(eigenGenes) <- retMods

    # subsetting modSz to match eigenGenes
    # For dropping communities within an iteration we use size and keep the larger since we have dynamic tree cut which uses topology.
    eigenGenes <- dropModuels(eigenGenes = eigenGenes,
                             Kurts = modSz[names(modSz) %in% retMods],
                              corCut = corCut)
    # when we drop between rounds we use kurtosis since we don't have access to topology at that point.
    tPc <- stats::prcomp(t(tEx))
    message(summary(tPc)$importance[2, 1:3])
    corPC <- stats::cor(tPc$x[, 1], t(eigenGenes))
    # order by cor with PC1 so that we regress out the eigengene most strongly associated with PC1.
    eigenGenes <- eigenGenes[order(abs(corPC), decreasing = TRUE), ]
  }
  return(eigenGenes)
}



####################################################################################################
# post network construction functions
####################################################################################################


#' compute_eigengene_matrix
#'
#' @param ex gene expression matrix with genes as rows and samples as columns
#' @param membership_matrix, a community membership (kME) matrix with genes as rows and communities as columns.
#' @param cutoff number of top genes to use when computing community signatures
#' @param pc_flag indicator. TRUE (default) means to use the 1st principal component (corrected for direction). FALSE uses the mean of scaled and centered top genes.
#'
#' @return A matrix with rows being the community signature and columns being samples
#'
#' @details Computes the community signatures (eigengenes) for an expression matrix given a particular community membership (kME) matrix. This
#' can be used to compute community signatures in a new expression dataset.
#' Note, community signatures are not corrected by icWGCNA iterations so it will not match signatures output from icWGCNA
#' if it is run on the network construction dataset.
#' When using these community signatures for modeling it may be best to include interaction terms or use tree based
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
#' @param t_memb `community_membership` or `full_community_membership` values from [icwgcna()]
#' @param K cutoff for top community genes to include for computing enrichment. Used in an AND condition with memb_cut.
#' @param memb_cut cutoff as a membership score threshold for determining top community genes for computing enrichment.  Used in an AND condition with K.
#' @param pangDB panglaoDB cell markers database. Default is to read data from the url [pangDB_link]
#' @param prolif list of proliferation genes to check. Default is to use [prolif_names]
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
                                         prolif = prolif_names)
{

  enr <- plyr::adply(t_memb,2, function(x) {
    c.types <- as.vector(stats::na.omit(unique(pangDB$`cell type`)))
    prolif_overlap <- table(rank(-x) <= K & x > memb_cut, rownames(t_memb) %in% prolif)
    ret <- t(plyr::ldply(c.types,function(ct) {
      if (ncol(prolif_overlap) > 1 &
          (prolif_overlap[2,2] > floor(.33 * length(prolif)))) {
        return(1)
      }

      overlap <- table(rank(-x) <= K & x > memb_cut,
                       rownames(t_memb) %in% pangDB$`official gene symbol`[pangDB$`cell type` == ct])
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


