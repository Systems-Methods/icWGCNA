#' icwgcna
#'
#' Iterative Correcting Weighted Gene Co-expression Network Analysis function constructing a network from an expression matrix.
#'
#' @param ex matrix of bulk RNA-seq or microarray gene expression data. This should be in log space and greater than 0.
#' @param expo exponent to use for soft thresholding. If NULL will use angular distance
#' @param Method correlation to use for distance measure, "pearson" (default) or "spearman"
#' @param q quantile (0-1) for first round filtering based on mean expression and standard deviation
#' @param maxIt maximum number of iterations must be 25 or less
#' @param maxComm maximum number of communities to be found
#' @param corCut correlation threshold used for dropping communities
#' @param covCut coefficient of variation (CoV) quantile threshold to use at each iteration  for selecting genes to build network. covCut = .667 would use the top third of genes based on CoV after regressing out largest community
#' @param mat_mult_method method for large matrix multiplication, "Rfast" (default) or "RcppEigen" (see {details})
#'
#' @return Returns a list with the following items:
#' * `community_membership` - community membership score (kME). Analogous to loadings in PCA.
#' * `community_signature` - community eigengene, the first principal component the expression of genes in this community. This can be thought of as the average of the scaled expression of top community genes.
#' * `full_community_membership` -  similar to community_membership but includes communities that were dropped in the iterative process.
#' * `full_community_signature` - similar to community_signature but includes communities that were dropped in the iterative process.
#' * `controlled_for` - The communities whose signatures were regressed out at each iteration.
#'
#' @details Iterative Correcting Weighted Gene Co-expression Network Analysis function for constructing a gene network from a gene expression matrix. The algorithm:
#'
#' 1. Constructs a signed wgcna network
#' 2. Drops correlated modules based on kurtosis.
#' 3. Regresses out the largest community from the expression data.
#' 4. Repeats steps 1-3 until a maximum number of communities or iterations is reached.
#'
#' Some differences from standard WGNCA (Horvath/Langfelder)
#'
#' - Makes heavy use of [Rfast][Rfast::Rfast-package] to compute adjacencies and TOM to enable iterative network creation on > 20K features.
#' - Uses signed adjacency in order to avoid possible distortions of community signatures (eigengenes).
#' - Iteratively regresses out strongest community in order to facilitate discovery of communities possibly obscured larger module(s).
#' - Clustering does not focus on merging communities but dropping to identify strongest module(s).
#' - Enables Spearman correlation for constructing adjacency matrix instead of Pearson to enable robust application in RNA-seq and micro-array data. Future updates may include mutual information
#'
#' For matrix multiplication the option "Rfast" will use [Rfast::mat.mult()],
#' which takes advantage of parallel processing across multiple cores. The option
#' "RcppEigen" will use the RcppEigen engine for C++ code, which tends to be faster
#' when using a single core, but does not take advantage of parallel processing across
#' multiple cores. If running this on a cluster with access to many computer core
#' there is a significant performance advantage to using [Rfast::mat.mult()]
#'
#'
#' @examples
#'
#'\dontrun{
#' library("UCSCXenaTools")
#' luad <- getTCGAdata(project = "LUAD", mRNASeq = TRUE, mRNASeqType = "normalized",
#'                   clinical = FALSE, download = TRUE)
#' ex <- data.table::fread(luad$destfiles)
#' results <- icwgcna(ex[,-1])}
#'
#' @export
icwgcna <- function(ex,
                    expo = 6,
                    Method = c("pearson", "spearman"),
                    q = .3,
                    maxIt = 10,
                    maxComm = 100,
                    corCut = .8,
                    covCut = .33,
                    mat_mult_method = c('Rfast', 'RcppEigen')) {
  # param checking
  if (!all(unlist(lapply(ex, is.numeric)))) {
    stop("all 'ex' columns must be numeric")
  }
  if (maxIt > 25 | maxIt < 1) {
    stop("maxIt must be between 1 and 25")
  }
  if (q >= 1 | q <= 0) {
    stop("q must be > 0 and <1 ")
  }
  if (corCut >= 1 | corCut <= 0) {
    stop("corCut must be > 0 and < 1")
  }
  if (covCut >= 1 | covCut <= 0) {
    stop("covCut must be > 0 and < 1")
  }

  Method <- match.arg(Method)
  mat_mult_method <- match.arg(mat_mult_method)


  SD <- matrix(NA, nrow(ex), maxIt)
  # standard deviation of each gene
  SD[, 1] <- apply(ex, 1, stats::sd)

  # need to remove any 0 SD genes
  SD_zero_index <- SD[, 1] == 0
  if (any(SD_zero_index)) {
    message('Removing ', sum(SD_zero_index),
            ' genes with a 0 standard deviation')
    ex <- ex[!SD_zero_index, ]
    SD <- SD[!SD_zero_index, ]
  }

  # average expression of each gene
  M <- apply(ex, 1, mean)
  # identify genes that should simply not be part of the first round due to low signal
  CoV      <- matrix(NA,nrow(ex), maxIt); CoV[,1] <- abs(SD[,1]/M)
  leaveOut <- M < stats::quantile(M, q) | SD[, 1] < stats::quantile(SD[, 1], q)
  tEx <- ex
  cont_for <- c()

  for (i in 1:maxIt) {
    # returns eigengenes for modules, orderd by |correlation to first pc of subsetted expression|
    tEigenGenes <- simpWGCNAsubNet(tEx[!leaveOut, ],
                                   expo = expo,
                                   Method = Method,
                                   n = 5,
                                   corCut = corCut,
                                   mat_mult_method = mat_mult_method
    )

    if (is.null(tEigenGenes)) {
      message("No more modules to be added. -- Stopping Iterations --")
      break()
    }

    rownames(tEigenGenes) <- paste0(LETTERS[i], 1:nrow(tEigenGenes))
    tMetaGenes <- as.data.frame(stats::cor(t(tEx), t(tEigenGenes), method = Method))

    if (i == 1) {
      eigenGenes <- tEigenGenes
      metaGenes <- tMetaGenes
      full_eigenGenes <- tEigenGenes
      full_metaGenes <- tMetaGenes
    } else {
      eigenGenes <- rbind(eigenGenes, tEigenGenes)
      metaGenes <- cbind(metaGenes, tMetaGenes)
      full_eigenGenes <- rbind(full_eigenGenes, tEigenGenes)
      full_metaGenes <- cbind(full_metaGenes, tMetaGenes)
      # not really size of community but more a measure of tightness
      clust_kurt <- apply(metaGenes, 2, Rfast::kurt)
      eigenGenes <- dropModuels(eigenGenes = eigenGenes,
                                Kurts = clust_kurt,
                                corCut = corCut,
                                verbose = FALSE)
      metaGenes <- metaGenes[, colnames(metaGenes) %in% rownames(eigenGenes)]
    }

    # this is the iterative correcting/controlled step were we regress out the largest signal (the eigengene from the largest module)
    if (i < maxIt) {
      cont_for <- c(cont_for, row.names(tEigenGenes)[1])
      # regress out 1st/largest eigen gene from this iteration
      x <- tEigenGenes[1, ]
      tEx <- plyr::aaply(.data = as.matrix(tEx),
                         .margins = 1,
                         .fun = function(y) {
                           stats::residuals.lm(stats::lm(y ~ x)) + mean(y)
                         })
      # filter on top third of genes based on coefficient of variation
      CoV[, i + 1] <- apply(tEx, 1, function(x) {
        abs(stats::sd(x) / mean(x))
      })
      # identify genes that should not be used to build the next sub-network due to low signal
      leaveOut <- CoV[, i] < stats::quantile(CoV[, i], covCut)
    }

    message(paste("Done with iteration:", i,
              ": current number of gene communities is",
              nrow(eigenGenes), "\n\n"))
    if (nrow(eigenGenes) >= maxComm) {
      message(paste("Number of modules", nrow(eigenGenes),
                  "is >=", maxComm,
                  "the maximum number of modules. -- Stopping Iterations --"))
      break()
    }
    if (i == maxIt) {
      message("Reached maximimum number of iterations")
      break()
    }
  }

  colnames(metaGenes) <- paste0("m", colnames(metaGenes))
  return(
    list(community_membership = metaGenes,
         community_signature = eigenGenes,
         full_community_membership = full_metaGenes,
         full_community_signature = full_eigenGenes,
         controlled_for = cont_for)
  )
}
