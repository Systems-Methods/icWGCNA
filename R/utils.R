
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
RfastTOMdist <- function(A,
                         mat_mult_method = c("Rfast", "RcppEigen")) {
  mat_mult_method <- match.arg(mat_mult_method)

  diag(A) <- 0
  A[is.na(A)] <- 0
  kk <- Rfast::colsums(A)
  denomHelp <- matrix(kk, ncol = length(kk), nrow = length(kk))
  denomTOM <- Rfast::Pmin(denomHelp, Rfast::transpose(denomHelp)) + (1 - A)
  rm(denomHelp, kk)

  if (mat_mult_method == "Rfast") {
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
fastTOMwrapper <- function(X,
                           expo = 6,
                           Method = c("pearson", "spearman"),
                           mat_mult_method = c("Rfast", "RcppEigen")) {
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
  modules <- dynamicTreeCut::cutreeHybrid(
    dendro = dend,
    distM = d,
    cutHeight = cutHeight,
    minClusterSize = 5,
    pamStage = TRUE,
    pamRespectsDendro = FALSE,
    verbose = 0
  )
  return(modules)
}

# drop correlated modules/communities based on kME kurtosis (keeping the community with higher kurtosis)
# Kurts is the metric used for selecting which community to drop
# it can be the size of each community, the average kME of each top X gene in each community or kurtosis each communities kME
# we use kurtosis
dropModuels <- function(eigenGenes,
                        Kurts = NULL,
                        corCut = .95,
                        method = c("pearson",
                                   "kendall",
                                   "spearman"),
                        logFlag=F) {
  method <- match.arg(method)
  if(logFlag){eigenGenes <- log(eigenGenes+1)}
  eigen_Cors <- stats::cor(t(eigenGenes), method = method)
  diag(eigen_Cors) <- 0
  while (any(eigen_Cors > corCut) & ncol(eigen_Cors) > 2) {
    doubleBreak <- FALSE

    for (i in 1:(nrow(eigen_Cors)))
    {
      for (j in (i + 1):nrow(eigen_Cors))
      {
        if (eigen_Cors[i, j] > corCut) {
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
          eigen_Cors <- stats::cor(t(eigenGenes), method = method)
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
  message(paste(
    "eigegenes trimmed to", nrow(eigenGenes),
    "due to correlation >", corCut,
    "max eigenCor =", max(signif(eigen_Cors, 2))
  ))
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
                            mat_mult_method = c("Rfast", "RcppEigen")) {
  Method <- match.arg(Method)
  mat_mult_method <- match.arg(mat_mult_method)

  message(paste(
    "Computing", nrow(tEx),
    "x", nrow(tEx),
    "TOM distance for subset of genes with higher variance"
  ))

  TOMd <- fastTOMwrapper(tEx,
    expo = expo,
    Method = Method,
    mat_mult_method = mat_mult_method
  )
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
    eigenGenes <- as.matrix(plyr::ldply(
      .data = retMods,
      .fun = function(m) {
        eg <- calcEigenGene(tEx[mods == m, ])
        return(eg)
      }
    ))
    rownames(eigenGenes) <- retMods

    # subsetting modSz to match eigenGenes
    # For dropping communities within an iteration we use size and keep the larger since we have dynamic tree cut which uses topology.
    eigenGenes <- dropModuels(
      eigenGenes = eigenGenes,
      Kurts = modSz[names(modSz) %in% retMods],
      corCut = corCut
    )
    # when we drop between rounds we use kurtosis since we don't have access to topology at that point.
    tPc <- stats::prcomp(t(tEx))
    message(summary(tPc)$importance[2, 1:3])
    corPC <- stats::cor(tPc$x[, 1], t(eigenGenes))
    # order by cor with PC1 so that we regress out the eigengene most strongly associated with PC1.
    eigenGenes <- eigenGenes[order(abs(corPC), decreasing = TRUE), ]
  }
  return(eigenGenes)
}
