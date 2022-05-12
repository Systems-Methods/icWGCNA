
####################################################################################################
# icWGCNA network construction related functions
####################################################################################################

# simple function computing angular distance to use as a adjacency measure.
angularDist <- function(x) {
  aDist <- asin(x) / (pi / 2)
  return(aDist)
}

# Wrapper of Rfast cora function to include spearman method which is used in calculating adjacency.
# x is a gene expression matrix with rows = genes and columns = samples.
RfastCor_wrapper <- function(x,
                             Method = "pearson") {
  gene_names <- rownames(x)
  x <- as.matrix(x)
  if (Method == "spearman") {
    x <- Rfast::rowRanks(x)
  }

  Cor <- Rfast::cora(t(x))

  rownames(Cor) <- gene_names
  colnames(Cor) <- gene_names

  return(Cor)
}


# TOM: topological overlap map from Ravasz, E., Somera, A., Mongru, D., Oltvai, Z. and Barab´asi, A. (2002). Science
# using Rfast package functions to  speed things up since we will be running wgcna up to 25 times
# A is a weighted adjacenty matrix
RfastTOMdist <- function(A) {
  diag(A) <- 0
  A[is.na(A)] <- 0
  kk <- Rfast::colsums(A)
  denomHelp <- matrix(kk, ncol = length(kk), nrow = length(kk))
  denomTOM <- Rfast::Pmin(denomHelp, Rfast::transpose(denomHelp)) + (1 - A)
  rm(denomHelp, kk)
  gc()

  # Rfast::mat.mult appear to be unstable in current release, coded in Rcpp from stack exchange
  # numTOM <- eigenMapMatMult(A,A) + A;
  numTOM <- Rfast::mat.mult(A, A) + A
  rm(A)
  gc()


  out <- 1 - as.matrix(numTOM / denomTOM)
  rm(denomTOM)
  gc()
  rm(numTOM)
  gc()
  diag(out) <- 0
  return(out)
}


# Topological overlap map (TOM) wrapper using Rfast package to speed up calculations,
# X is expression matrix w each column being one sample
# This funcion can compute TOM for 24K gene matrix in 8 minute of an AWS-EC2 c5.18xlarge instance,
# though in practice we run it on subsets of most varialbe genes
# note that we only use signed and weighted adjacencies
fastTOMwrapper <- function(X, expo = NULL, Method = "pearson") {
  # compute weighted adjacency matrix using Rfast package
  X <- RfastCor_wrapper(X,
                        Method = Method)
  X[is.na(X)] <- 0
  X[X < 0] <- 0

  if (is.null(expo)) {
    A <- angularDist(X)
  } else {
    A <- X^expo
    A[X < 0] <- 0
  }
  rm(X)
  gc()

  tom <- RfastTOMdist(A)
  rm(A)
  gc()
  return(tom)
}


# computes the eigengene for a single community
calcEigenGene <- function(tEx) {
  if (nrow(tEx) == 1) {
    return(unlist(tEx))
  }
  tEx <- as.matrix(tEx)
  pc1 <- prcomp(t(tEx))$x[, 1]
  if (sum(cor(t(tEx), pc1) < 0, na.rm = TRUE) / nrow(tEx) > .5) {
    pc1 <- -pc1
  }
  return(pc1)
}


# wrapper of cutreeHybrid from Langfelder and Horvath's dynamic tree cutting pacakge.(https://cran.r-project.org/web/packages/dynamicTreeCut/)
# Instead of merging modules like they do, we'll be dropping correlated modules based on kME kurtosis
# i.e. we'll be selecting between two correlated modules
cutreeHybridWrapper <- function(d,
                                quantCut = 0.75) {
  dend <- hclust(as.dist(d), method = "average")
  refHeight <- quantile(dend$height, .05, type = 1)
  cutHeight <- as.numeric(quantCut * (max(dend$height) - refHeight) + refHeight)
  modules <- dynamicTreeCut::cutreeHybrid(dendro = dend,
                                          distM = d,
                                          cutHeight = cutHeight,
                                          minClusterSize = 5,
                                          pamStage = TRUE,
                                          pamRespectsDendro = FALSE)
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
  eigen_Cors <- cor(t(eigenGenes))
  diag(eigen_Cors) <- 0
  while (any(eigen_Cors > corCut)) {
    doubleBreak <- FALSE

    for (i in 1:(nrow(eigen_Cors)))
    {
      for (j in (i + 1):nrow(eigen_Cors))
      {
        if (eigen_Cors[i, j] > corCut) {
          if (verbose) {
            print(paste("threshold", corCut, "exceeded for mods",
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
          eigen_Cors <- cor(t(eigenGenes))
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
  eigen_Cors <- cor(t(eigenGenes))
  diag(eigen_Cors) <- 0
  print(paste("eigegenes trimmed to", nrow(eigenGenes),
              "due to correlation >", corCut,
              "max eigenCor =", max(signif(eigen_Cors, 2))))
  return(eigenGenes)
}


# This is the simplified, Rfast-based WGCNA function that is used in each iteration
# Rfast is used to speed up both the adjacency calculation and the TOM computation
simpWGCNAsubNet <- function(tEx,
                            expo = 6,
                            Method = "pearson",
                            n = 15,
                            minMods = 5,
                            corCut = .6) {
  print(paste("Computing", nrow(tEx),
              "x", nrow(tEx),
              "TOM distance for subset of genes with higher variance"))

  TOMd <- fastTOMwrapper(tEx,
                         expo = expo,
                         Method = Method)
  gc()
  mods <- cutreeHybridWrapper(TOMd)$labels
  modSz <- table(mods)
  if (sum(modSz >= n) < minMods) {
    retMods <- names(modSz[rank(-modSz) <= minMods])
  } else {
    retMods <- names(modSz[modSz >= n])
  }

  # omit the catch all module for genes that do not belong to a community
  retMods <- retMods[retMods != "0"]
  print(paste("number of modules found is", length(retMods)))

  eigenGenes <- as.matrix(plyr::ldply(.data = retMods,
                                      .fun = function(m) {
                                        eg <- calcEigenGene(tEx[mods == m, ])
                                        return(eg)
                                      }))
  rownames(eigenGenes) <- retMods

  # subsetting modSz to match eigenGenes
  eigenGenes <- dropModuels(eigenGenes = eigenGenes,
                            Kurts = modSz[names(modSz) %in% retMods],
                            corCut = corCut)
  tPc <- prcomp(t(tEx))
  print(summary(tPc)$importance[2, 1:3])
  corPC <- cor(tPc$x[, 1], t(eigenGenes))
  eigenGenes <- eigenGenes[order(abs(corPC), decreasing = TRUE), ]

  gc()
  return(eigenGenes)
}



####################################################################################################
# post network construction functions
####################################################################################################


#' compute_eigengene_matrix
#'
#' @param ex gene expression matrix with genes as rows and samples as columns
#' @param membership_matrix, a community membership (kME) matrix with genes as rows and commutities as columns.
#' @param cutoff number of top genes to use when computing community signatures
#' @param pc_flag indicator. T (default) means to use the 1st principal component (corrected for direction). FALSE uses the mean of scaled and centered top genes.
#'
#' @return A matrix with rows being the community signature and columns being samples
#' @export
#'
#' @details Computes the community signatures (eigengenes) for an expression matrix given a particular community membership (kME) matrix. This
#' can be used to compute community signatures in a new expression dataset.
#' Note, community signatures are not corrected by icWGCNA iterations so it will not match signatures output from icWGCNA
#' if it is run on the network constuction dataset.
#' When using these community signatures for modeling it may be best to include interaction terms or use tree based
#' methods since dependencies are not addressed in this output matrix.
#'
#' @examples
compute_eigengene_matrix <- function(ex,
                                     membership_matrix,
                                     cutoff = 5,
                                     pc_flag = TRUE) {
  ex <- ex[apply(as.matrix(ex), 1, sd) != 0, ]
  meta <- meta[rownames(meta) %in% rownames(ex), ]
  m_genes <- rownames(meta)
  e_genes <- rownames(ex)

  tEigen <- plyr::aaply(as.matrix(meta), 2, function(k) {
    tInds <- rank(-k) <= cutoff
    tEx <- ex[e_genes %in% m_genes[tInds], ]
    if (all(tEx == 0)) {
      return(rep(0, ncol(tEx)))
    }

    scaled_ave <- apply(t(scale(t(tEx))), 2, mean)
    if (pc_flag == FALSE) {
      return(scaled_ave)
    } else {
      pc1 <- prcomp(t(tEx))$x[, 1]
      if (cor(scaled_ave, pc1) < 0) {
        pc1 <- -pc1
      }
      return(pc1)
    }
  })
}
