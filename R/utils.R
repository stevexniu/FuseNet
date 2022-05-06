#' @useDynLib FuseNet
#' @importFrom Rcpp sourceCpp
#' 
NULL

#' Data Normalization
#'
#' Normalize count or other data.
#'
#' @param counts Data to normalize. An N x M matrix with N rows of features and M columns of data points.
#' @param normalization Normalization method used. Default is cosine.
#' \itemize{
#'   \item cosine, cosine normalization: feature counts for each data point are divided by the L2 norm of them.
#'   \item lognorm, log normalization: feature counts for each data point are divided by the total sum of them. Then the data is multiplied by the scale.factor before taking a log-transformed by log(1+x).
#'   \item none, additional normalization is not performed.
#'   }
#' @param normalize_factor Normalization factor used with lognorm method. Default is 10000.
#' @param zero_percent Zero-entry percentage threshold. If the number of zeros in the returned matrices is above this number, a sparse matrix will be returned. Default is 0.7 aka 70\%.
#' @param verbose Whether to display a process bar. Default is FALSE.
#' @return Returns the normalized data.
#' @importFrom Matrix nnzero Matrix
#' @examples \dontrun{
#' normalized.data <- Normalization(data)
#' }
Normalization <- function(counts, normalization = c("cosine", "lognorm", "none"), normalize_factor = 1e4, zero_percent = 0.7, verbose = FALSE){
  if(nnzero(counts)/length(counts) < (1-zero_percent)){
    counts <- Matrix(data = counts, sparse = TRUE)
    normalization <- match.arg(arg = normalization)
    return(switch(normalization,
                  lognorm = LogNormSparse(data = counts, scale_factor = normalize_factor, display_progress = verbose),
                  cosine = CosineNormSparse(data = counts, display_progress = verbose),
                  none = counts))
  } else {
    counts <- as.matrix(x = counts)
    normalization <- match.arg(arg = normalization)
    return(switch(normalization,
                  lognorm = LogNorm(data = counts, scale_factor = normalize_factor, display_progress = verbose),
                  cosine = CosineNorm(data = counts, display_progress = verbose),
                  none = counts))
  }
}

#' Data Scaling
#'
#' Scale and center data.
#'
#' @param matrix Matrix to use.
#' @param verbose Whether to display a process bar. Default is FALSE.
#' @return Returns a scaled and centered data.
#' @examples \dontrun{
#' scaled.data <- Scaling(data)
#' }
Scaling <- function(matrix, verbose = FALSE){
  if(is(matrix,"matrix")){
    scale.data <- FastRowScale(mat = matrix, display_progress = verbose)
  } else if (is(matrix,"dgCMatrix")){
    scale.data <- FastSparseRowScale(mat = matrix, display_progress = verbose)
  }
  dimnames(scale.data) <- dimnames(matrix)
  return(scale.data)
}

#' Principal Components Analysis
#'
#' Perform principal components analysis on scaled and centered data (z-scores).
#'
#' @param scaled_data Scaled and centered data (z-scores).
#' @param n_dims Number of dimensions to return. Default is 10.
#' @param seed Random seed number. Default is 1.
#' @return Returns a list contains outputs from \code{\link[irlba]{irlba}}.
#' @importFrom irlba irlba
#' @examples \dontrun{
#' pca.result <- PCA(scaled.data, n_dims = 3)
#' }
PCA <- function(scaled_data, n_dims = 10, seed = 1){
  set.seed(seed = seed)
  pca.res <- irlba(A = scaled_data, nv = n_dims)
  return(pca.res)
}

#' L1-like Norm
#'
#' Compute L1-like norm defined by the absolute values of the differences between each entry of two matrices with the same dimensions.
#'
#' @param mat1,mat2 Two matrices.
#' @return Returns a matrix with L1-like norms.
#' @examples \dontrun{
#' l1.norm <- L1Norm(matrix1, matrix2)
#' }
L1Norm <- function(mat1, mat2){
  return(abs(x = (mat1 - mat2)))
}

#' L2-like Norm
#'
#' Compute L2-like norm defined by the square values of the differences between each entrie of two matrices with the same dimensions.
#'
#' @param mat1,mat2 Two matrices.
#' @return Returns a matrix with L2-like norms.
#' @examples \dontrun{
#' l2.norm <- L2Norm(matrix1, matrix2)
#' }
L2Norm <- function(mat1, mat2){
  return((mat1 - mat2)^2)
}

#' Euclidean Distance
#'
#' Compute euclidean nearest neighbor distances.
#'
#' @param data An M x d matrix with M rows of data points and d columns of features.
#' @param ka Number of nearest neighbors. See details from \code{\link[RANN]{nn2}}.
#' @return Returns the distance matrix.
#' @importFrom RANN nn2
#' @examples \dontrun{
#' dist.mat <- EuclideanDist(data, ka = 10)
#' }
EuclideanDist <- function(data, ka){
  knn <- nn2(data = data, k = ka)
  n <- nrow(x = knn[[1]])
  dist.mat <- matrix(data = 0, nrow = n, ncol = n)
  for(i in 1:nrow(x = knn[[1]])){
    dist.mat[i, knn[[1]][i, ]] <- knn[[2]][i, ]
  }
  return(dist.mat)
}

#' Gaussian Distance
#'
#' Compute gaussian nearest neighbor distances.
#'
#' @param data An M x d matrix or data.frame with M rows of data points and d columns of features.
#' @param ka Number of nearest neighbors. See details from \code{\link[RANN]{nn2}}.
#' @return Returns the distance matrix.
#' @importFrom RANN nn2
#' @examples \dontrun{
#' dist.mat <- GaussianDist(data, ka = 10)
#' }
GaussianDist <- function(data, ka){
  knn <- nn2(data = data, k = ka * 3)
  n <- nrow(x = knn[[1]])
  dist.mat <- matrix(data = 0, nrow = n, ncol = n)
  for(i in 1:nrow(x = knn[[1]])) {
    dist.mat[i, knn[[1]][i, ]] <- exp(x = (-1 * (knn[[2]][i, ] / knn[[2]][i, ka])^2))
  }
  dist.mat <- dist.mat + t(x = dist.mat)
  dist.mat <- sweep(x = dist.mat, MARGIN = 1, STATS = rowSums(x = dist.mat), FUN = "/")
  dist.mat[is.na(x = dist.mat)] <- 0
  return(dist.mat)
}

#' Run Geometric Sketching
#'
#' Run Geometric sketching sampling.
#' See \href{https://www.sciencedirect.com/science/article/pii/S2405471219301528}{Hie et al., 2019} for details.
#' This function relies on reticulate R package to import geosketch python library. Only one python environment is allowed by reticulate configuration.
#' Therefore please be careful when importing other R packages that will invoke and set up reticulate python configuration, i.e library(Seurat), before calling GeomSketch function.
#'
#' @param data An M x d matrix or data.frame with M rows of data points and d columns of features.
#' @param geom_size Size of geometric sketches to return.
#' @param is_pca Whether the data columns are principal components.
#' @param n_pca Number of PCA dimensions to use. Default is 10.
#' @param seed Random seed number. Default is 1.
#' @param which_python Path to python3 used.
#' @import reticulate
#' @importFrom irlba irlba
#' @return Returns geometric sketch ID.
#' @export
#' @examples \dontrun{
#' sketch.ids <- RunGeomSketch(data, geom_size = 1000)
#' }
RunGeomSketch <- function(data, geom_size, is_pca = FALSE, n_pca = 10, seed = 1, which_python = Sys.which(names = "python3")){
  use_python(python = which_python, required = TRUE)
  if(!py_module_available("geosketch")){
    stop("Cannot find geosketch, please install through pip (e.g. pip install geosketch).")
  }
  geosketch <- import(module = 'geosketch')
  if(is_pca){
    sketch.ind <- unlist(x = geosketch$gs(data$u, as.integer(x = geom_size))) + 1
  } else {
    set.seed(seed = seed)
    pca.res <- irlba(A = data, nv = n_pca)
    sketch.ind <- unlist(x = geosketch$gs(pca.res$u, as.integer(x = geom_size))) + 1
  }
  return(sketch.ind)
}

#' Bootstrapping of Distance Matrix
#'
#' Bootstrapping  of Distance Matrix
#'
#' @param data An M x d matrix or data.frame with M rows of data points and d columns of features.
#' @param dist_mat_null An M x M distance matrix calculated from the original data (null).
#' @param k Number of nearest neighbors. Default is 10. See details from \code{\link[RANN]{nn2}}.
#' @param kernel Kernel distance used:
#' \itemize{
#'   \item gaussian, gaussian distance kernel. See details \code{\link[FuseNet]{EuclideanDist}}.
#'   \item euclidean, euclidean distance kernel. See details \code{\link[FuseNet]{GaussianDist}}.
#'   }
#' @param normalization Normalization method used:
#' #' \itemize{
#'   \item cosine, cosine normalization. See details \code{\link[FuseNet]{Normalization}}.
#'   \item lognorm, log normalization. See details \code{\link[FuseNet]{Normalization}}.
#'   \item none, normalization is not performed.
#'   }
#' @param normalize_factor Normalize factor used in log normalization. Default is 10000. See details \code{\link[FuseNet]{Normalization}}.
#' @param pca_dims Number of dimensions used. Default is 0 and PCA is not performed.
#' @param norm_type Type of norm used:
#' \itemize{
#'   \item l1, L1-like norm. See details \code{\link[FuseNet]{L1Norm}}.
#'   \item l2, L1-like norm. See details \code{\link[FuseNet]{L2Norm}}.
#'   }
#' @param n_iters Number of bootstrapping iterations. Default is 100.
#' @param ratio Fraction of features to be downsampled in the original data matrix. Default is 0.05 aka 5\%.
#' @param t Matrix power used for the distance matrix. Default is 0 and powering is not performed. See \code{\link[FuseNet]{MatrixPower}} for details.
#' @param calc_perturb_mat Whether to calculate the perturb matrix. Default is FALSE.
#' @param n_cores Number of cores used. Default is to use all existing cores. See details \code{\link[parallel]{makeCluster}}.
#' @param zero_percent Zero-entry percentage threshold. If the number of zeros in the returned matrices is above this number, a sparse matrix will be returned. Default is 0.7 aka 70\%.
#' @param ... Additional parameters pass to \code{\link[parallel]{makeCluster}}.
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom Matrix nnzero Matrix
#' @import foreach
#' @return Returns a list with entries:
#' \itemize{
#'   \item feature_weight, n x d binary matrix with n rows of bootstrap iterations and d columns of features where 0 means feature not sampled and 1 means sampled.
#'   \item sample_weight, n x M matrix with n rows of bootstrap iterations and M columns of data points where each entry represents weight.
#'   \item perturb_mat, d x M matrix with d rows of features and M columns of data points where each entry represents the relative importance of a feature to a data point.
#'   \item dist_mat, M x M distance matrix.
#'   }
#' @examples \dontrun{
#' boot.result <- Bootstrap(data, dist.mat.null, k = 10, pca_dims = 10, n_iters = 100, ratio = 0.1)
#' }
Bootstrap <- function(data, dist_mat_null, k = 10, kernel = c("gaussian", "euclidean"), normalization = c("cosine", "lognorm", "none"), normalize_factor = 1e4, pca_dims = 0, norm_type = c("l1", "l2"), 
                      n_iters = 100, ratio = 0.05, t = 0, calc_perturb_mat = FALSE, n_cores = NULL, zero_percent = 0.7, ...){
  kernel <- match.arg(arg = kernel)
  dist.fxn <- switch(kernel,
                     gaussian = GaussianDist,
                     euclidean = EuclideanDist)
  norm_type <- match.arg(arg = norm_type)
  norm.fxn <- switch(norm_type,
                     l1 = L1Norm,
                     l2 = L2Norm)
  n_cell <- nrow(x = data)
  n_feature <- ncol(x = data)
  n_sample <- n_feature * ratio
  dist_mat_null <- as.matrix(dist_mat_null)
  
  if(is.null(x = n_cores)) n_cores <- detectCores() - 1
  cl <- makeCluster(spec = n_cores, ...)
  registerDoParallel(cl = cl)
  if(pca_dims > 0){
    results <- foreach(i = 1:n_iters, .multicombine = TRUE, .packages = "FuseNet") %dopar% {
      set.seed(seed = i)
      sample.idx <- sample(x = 1:n_feature, size = n_sample)
      data.temp <- PCA(scaled_data = data[, sample.idx], n_dims = pca_dims)
      dist.mat.temp <- dist.fxn(data = data.temp$u, ka = k)
      dist.mat.temp <- MatrixPower(A = as.matrix(x = dist.mat.temp), n = t)
      norm.temp <- norm.fxn(mat1 = dist.mat.temp, mat2 = dist_mat_null)
      list(sample.idx, colSums(x = norm.temp), norm.temp)
    }
  } else {
    results <- foreach(i = 1:n_iters, .multicombine = TRUE, .packages = "FuseNet") %dopar% {
      set.seed(seed = i)
      sample.idx <- sample(x = 1:n_feature, size = n_sample)
      dist.mat.temp <- dist.fxn(data = data[, sample.idx], ka = k)
      dist.mat.temp <- MatrixPower(A = as.matrix(x = dist.mat.temp), n = t)
      norm.temp <- norm.fxn(mat1 = dist.mat.temp, mat2 = dist_mat_null)
      list(sample.idx, colSums(x = norm.temp), norm.temp)
    }
  }
  stopCluster(cl = cl)
  gc()

  message("Finalize")
  feature.weight.mat <- lapply(X = results, FUN = "[[", ... = 1L)
  sample.weight <- lapply(X = results, FUN = "[[", ... = 2L)
  dist.mat <- lapply(X = results, FUN = "[[", ... = 3L)
  
  feature.weight.mat <- do.call(what = "rbind", args = feature.weight.mat)
  sample.weight <- do.call(what = "rbind", args = sample.weight)
  dist.mat <- Reduce(f = "+", x = dist.mat)

  feature.weight <- matrix(data = 0, nrow = n_iters, ncol = n_feature)
  for(i in 1:n_iters) {
    feature.weight[i, feature.weight.mat[i, ]] <- 1
  }

  if(nnzero(x = feature.weight)/length(x = feature.weight) < (1-zero_percent)){
    feature.weight <- Matrix(data = feature.weight, sparse = TRUE)
  } else {
    feature.weight <- as.matrix(x = feature.weight)
  }
  if(nnzero(x = sample.weight)/length(x = sample.weight) < (1-zero_percent)){
    sample.weight <- Matrix(data = sample.weight, sparse = TRUE)
  } else {
    sample.weight <- as.matrix(x = sample.weight)
  }
  if(nnzero(x = dist.mat)/length(x = dist.mat) < (1-zero_percent)){
    dist.mat <- Matrix(data = dist.mat, sparse = TRUE)
  } else {
    dist.mat <- as.matrix(x = dist.mat)
  }
  dimnames(x = feature.weight) <- list(paste("bootstrap", 1:n_iters, sep = "_"), colnames(x = data))
  dimnames(x = sample.weight) <- list(paste("bootstrap", 1:n_iters, sep = "_"), rownames(x = data))
  dimnames(x = dist.mat) <- dimnames(x = dist_mat_null)
  if(calc_perturb_mat){
    perturb.mat <- t(x = feature.weight) %*% sample.weight
  } else {
    perturb.mat <- matrix(nrow = 0,ncol = 0)
  }
  dist.mat <- (dist.mat + t(x = dist.mat)) / n_iters  
  return(list(feature_weight = feature.weight, sample_weight = sample.weight, perturb_mat = perturb.mat, dist_mat = dist.mat))
}

#' Project Geometric Sketches
#'
#' Project weights learned from geometric sketches onto other data points based on nearest neighbors.
#'
#' @param data An M x d matrix or data.frame with M rows of data points and d columns of features.
#' @param data_query An N x d matrix or data.frame with N rows of queried data points and d columns of features.
#' @param weights_data A vector of weights for each data point with total length of M.
#' @param k Number of nearest neighbors. See details from \code{\link[RANN]{nn2}}.
#' @importFrom RANN nn2
#' @return Returns the weight matrix of the queried data.
#' @examples \dontrun{
#' weight.mat <- ProjectSketch(data, data.query, weights.data, 30)
#' }
ProjectSketch <- function(data, data_query, weights_data, k){
  knn.obj <- nn2(data = data, query = data_query, k = k)
  gauss <- function(x, k) exp(x = (-1 * (x / x[k])^2))
  k.gauss <- ceiling(x = k / 3)
  knn.obj[[2]] <- t(x = apply(X = knn.obj[[2]], MARGIN = 1, FUN = function(x) gauss(x, k.gauss) / sum(gauss(x, k.gauss))))
  n.query <- nrow(x = data_query)
  n.data <- nrow(x = data)
  dist.mat <- matrix(data = 0, nrow = n.query, ncol = n.data)
  for(i in 1:n.query){
    dist.mat[i, knn.obj[[1]][i, ]] <- knn.obj[[2]][i, ]
  }
  weights.query <- dist.mat %*% weights_data
  return(weights.query)
}

#' Matrix Big Power
#'
#' Compute big powers (>10) of matrix.
#'
#' @param A A matrix.
#' @param n Number of powers.
#' @return Returns the powered matrix.
#' @examples \dontrun{
#' A15 <- MatrixBigPower(A, 15)
#' }
MatrixBigPower <- function(A, n){
  if(n > 0){
    e <- eigen(x = A)
    M <- e$vectors
    d <- e$values
    return(M %*% diag(x = d^n) %*% solve(a = M))
  } else {
    return(A)
  }
}

#' Matrix Power
#'
#' Compute powers of matrix.
#'
#' @param A A matrix.
#' @param n Number of powers.
#' @return Returns the powered matrix.
#' @importFrom expm expm %^%
#' @examples \dontrun{
#' A.cubic <- MatrixPower(A, 3)
#' }
MatrixPower <- function(A, n){
  if(n > 0){
    return(A %^% n)
  } else {
    return(A)
  }
}

#' Merge Two List
#'
#' Merge twp lists by the name of their elements.
#'
#' @param a,b Two lists.
#' @return Returns the merged list.
MergeLists <- function(a, b) {
  a.names <- names(a)
  b.names <- names(b)
  m.names <- unique(c(a.names, b.names))
  sapply(m.names, function(i) {
    if (i %in% b.names) b[[i]]
    else a[[i]]
  }, simplify = FALSE)
}
