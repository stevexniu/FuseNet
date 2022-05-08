#' @include utils.R
#'
NULL

#' Initialize FuseNet Object
#'
#' Initialize a FuseNet object.
#'
#' @param raw_data Raw data. An N x M matrix with N rows of features and M columns of data points.
#' @param project_name Project name. Default is none.
#' @param normalization Normalization method used. Default is cosine. See details \code{\link[FuseNet]{Normalization}}.
#' \itemize{
#'   \item cosine, cosine normalization: feature counts for each data point are divided by the L2 norm of them.
#'   \item lognorm, log normalization: feature counts for each data point are divided by the total sum of them. Then the data is multiplied by the scale.factor before taking a log-transformed by log(1+x).
#'   \item none, additional normalization is not performed.
#'   }
#' @param normalize_factor Normalization factor used with lognorm method. Default is 10000.
#' @param zero_percent Zero-entry percentage threshold. If the number of zero entries in the returned matrices is above this number, a sparse matrix will be returned. Default is 0.7 aka 70\%.
#' @param pca_dims Number of dimensions used. Default is 0 and PCA is not performed.
#' @param kernel Kernel distance used:
#' \itemize{
#'   \item gaussian, gaussian distance kernel. See details \code{\link[FuseNet]{EuclideanDist}}.
#'   \item euclidean, euclidean distance kernel. See details \code{\link[FuseNet]{GaussianDist}}.
#'   }
#' @param k Number of nearest neighbors. Default is 100. See details from \code{\link[RANN]{nn2}}.
#' @param t Matrix power used for the distance matrix. Default is 0 and powering is not performed. See \code{\link[FuseNet]{MatrixPower}} for details.
#' @param verbose Whether to display a process bar. Default is FALSE.
#' @param seed Random seed number. Default is 1.
#' @return Returns a FuseNet object.
#' @importFrom Matrix nnzero Matrix t
#' @importFrom methods is new
#' @importFrom stats sd
#' @export
#' @concept fuse
#' @examples {
#' object <- InitiateFuseNet(t(iris[,1:4]), project_name = "FuseNet", k = 10)
#' }
#' 
InitiateFuseNet <- function(raw_data, project_name = "", normalization = c("cosine", "lognorm", "none"), normalize_factor = 1e4, zero_percent = 0.7, pca_dims = 0, kernel = c("gaussian", "euclidean"), k = 100, t = 0, verbose = FALSE, seed = 1){
  message("Initiate FuseNet")
  feature.sd <- apply(X = raw_data, MARGIN = 1, FUN = sd)
  if(any(feature.sd == 0)){
    message("Removing Missing Features")
    raw_data <- raw_data[-which(x = feature.sd == 0),]
  }  
  object <- new(Class = "FuseNet", raw_data = raw_data, project_name = project_name)
  normalization <- match.arg(arg = normalization)
  if(normalization != "none") message("Normalize Data")
  object@normalized_data <- Normalization(counts = object@raw_data, normalization = normalization, normalize_factor = normalize_factor, zero_percent = zero_percent, verbose = verbose)
  dimnames(x = object@normalized_data) <- dimnames(x = object@raw_data)
  if(pca_dims > 0){
    message("Scale Data")
    object@scaled_data <- Scaling(matrix = object@normalized_data, verbose = verbose)
    dimnames(x = object@scaled_data) <- dimnames(x = object@raw_data)
    message("Run PCA")
    object@pca <- PCA(scaled_data = t(x = object@scaled_data), n_dims = pca_dims, seed = seed)
    data.use <- object@pca$u
  } else {
    data.use <- t(x = object@normalized_data)
  }
  message("Find Nearest Neighbors")
  kernel <- match.arg(arg = kernel)
  dist.mat <- switch(kernel,
                             euclidean = EuclideanDist(data = data.use, ka = k),
                             gaussian = GaussianDist(data = data.use, ka = k))
  message("Matrix Power")
  dist.mat <- MatrixPower(A = as.matrix(x = dist.mat), n = t)
  message("Finalize")
  if(nnzero(x = dist.mat)/length(x = dist.mat) < (1 - zero_percent)) dist.mat <- Matrix(data = dist.mat, sparse = TRUE)
  dimnames(x = dist.mat) <- list(colnames(x = object@raw_data), colnames(x = object@raw_data))
  object@dist_null <- dist.mat
  object@params <- mget(x = names(x = formals()), envir = sys.frame(which = sys.nframe()))[-1]
  return(object)
}

#' Geometric Sketching
#'
#' Run Geometric sketching sampling.
#' See Hie et al. (\doi{10.1016/j.cels.2019.05.003}) for details.
#'
#' @param object A FuseNet object.
#' @param geom_size Size of geometric sketches to return. Default is 1000.
#' @param geom_pca_dims Number of PCA dimensions to use. Default is 10.
#' @param sketch_k Number of nearest neighbors for the sketched data points. Default is 30. See details from \code{\link[FuseNet]{InitiateFuseNet}}.
#' @param sketch_t Matrix power used for the distance matrix. Default is 0 and powering is not performed. See details from \code{\link[FuseNet]{InitiateFuseNet}}.
#' @param sketch_n_pca Number of dimensions used for the sketched data points. Default is 0 and PCA is not performed. See details from \code{\link[FuseNet]{InitiateFuseNet}}.
#' @param geom_pca_seed Random seed number for PCA. Default is 1.
#' @param which_python Path to python3 used.
#' @return Returns a FuseNet object contains geometric sketch IDs.
#' @export
#' @concept fuse
#' @examples \dontrun{
#' object <- GeomSketch(object, geom_size = 10)
#' }
#' 
GeomSketch <- function(object, geom_size = 1000, geom_pca_dims = 10, sketch_n_pca = 0, sketch_k = 30, sketch_t = 0, geom_pca_seed = 1, which_python = Sys.which(names = "python3")){
  if(nrow(x = object@scaled_data) == 0){
    message("Scale Data")
    object@scaled_data <- Scaling(matrix = object@normalized_data, verbose = object@params$verbose)
  } 
  object@sketch_id <- RunGeomSketch(data = t(object@scaled_data), geom_size = as.integer(geom_size), is_pca = FALSE, n_pca = geom_pca_dims, seed = geom_pca_seed, which_python = which_python)
  params.use <- object@params
  object.temp <- InitiateFuseNet(object@raw_data[, object@sketch_id], project_name = params.use$project_name, normalization = params.use$normalization, normalize_factor = params.use$normalize_factor, zero_percent = params.use$zero_percent, pca_dims = sketch_n_pca, kernel = params.use$kernel, k = sketch_k, t = sketch_t, verbose = params.use$verbose, seed = params.use$seed)
  object@sketch_dist <- object.temp@dist_null
  object@params <- MergeLists(object@params, mget(x = names(x = formals()), envir = sys.frame(which = sys.nframe()))[-1])
  return(object)
}

#' Run Data Fusion
#'
#' Run Data Fusion.
#'
#' @param object A FuseNet object.
#' @param n_iters Number of bootstrapping iterations. Default is 100.
#' @param ratio Fraction of features to be downsampled in the original data matrix. Default is 0.05 aka 5\%.
#' @param pca_dims Number of principle components. Default is 0 and PCA is not run.
#' @param k Number of nearest neighbors used. Default is 100.
#' @param t Matrix power used for the distance matrix. Default is 0 and powering is not performed.
#' @param norm_type Type of norm used:
#' \itemize{
#'   \item l1, L1-like norm. See details \code{\link[FuseNet]{L1Norm}}.
#'   \item l2, L1-like norm. See details \code{\link[FuseNet]{L2Norm}}.
#'   }
#' @param return_perturb_mat Whether to return the perturb matrix. Default is FALSE.
#' @param n_cores Number of cores used. Default is to use all existing cores. See details \code{\link[parallel]{makeCluster}}.
#' @param ... Additional parameters pass to \code{\link[parallel]{makeCluster}}.
#' @return Returns a FuseNet object.
#' @export
#' @concept fuse
#' @examples {
#' object <- InitiateFuseNet(t(iris[,1:4]), project_name = "FuseNet", k = 3)
#' object <- RunFuseNet(object, n_iters = 1, k = 10, ratio = 0.5, n_cores = 1)
#' }
#' 
RunFuseNet <- function(object, n_iters = 100, ratio = 0.05, pca_dims = 0, k = 100, t = 0, norm_type = c("l1", "l2"), return_perturb_mat = FALSE, n_cores = NULL, ...){
  message("Run Bootstrapping")
  params.use <- object@params
  if(params.use$pca_dims > 0){
    data.use <- t(x = object@scaled_data)
  } else {
    data.use <- t(x = object@normalized_data)
  }
  if(length(x = object@sketch_id) > 0){
    data.use <- data.use[object@sketch_id, ]
    dist.use <- object@sketch_dist
  } else {
    dist.use <- object@dist_null
  }
  bootstap.results <- Bootstrap(data = data.use, dist_mat_null = dist.use, k = k, kernel = params.use$kernel, pca_dims = pca_dims, norm_type = norm_type, n_iters = n_iters, ratio = ratio, t = t, calc_perturb_mat = return_perturb_mat, n_cores = n_cores, zero_percent = params.use$zero_percent, ...)
  object@weight_mat <- bootstap.results[1:3]
  object@dist_mat <- bootstap.results$dist_mat
  object@params <- MergeLists(object@params, mget(x = names(x = formals()), envir = sys.frame(which = sys.nframe()))[-1])
  return(object)
}

#' Weight and Fuse Objects
#'
#' Weight and fuse distance matrices based on the relative weights.
#'
#' @param ... FuseNet objects to fuse.
#' @param project_k Number of nearest neighbors to project based on geometric sketches, if it has been run. Default is 10.
#' @param zero_percent Zero-entry percentage threshold. If the number of zero entries in the returned matrices is above this number, a sparse matrix will be returned. Default is 0.7 aka 70\%.
#' @return Returns a list with entries:
#' \itemize{
#'   \item fused_weight, n x M matrix with n rows of objects or modality to fuse and M columns of data points.
#'   \item fused_dist, M x M fused distance matrix.
#'   }
#' @importFrom Matrix nnzero Matrix t
#' @export
#' @concept fuse
#' @examples {
#' object1 <- InitiateFuseNet(t(iris[,1:2]), project_name = "FuseNet", k = 3)
#' object1 <- RunFuseNet(object1, n_iters = 1, k = 3, ratio = 0.5, n_cores = 1)
#' object2 <- InitiateFuseNet(t(iris[,3:4]), project_name = "FuseNet", k = 3)
#' object2 <- RunFuseNet(object2, n_iters = 2, k = 3, ratio = 0.5, n_cores = 1)
#' fused.data <- FuseData(object1, object2)
#' }
#' 
FuseData <- function(..., project_k = 10, zero_percent = 0.7){
  obj.list <- list(...)
  n.obj <- length(x = obj.list)
  n.cells <- unlist(lapply(X = obj.list, FUN = function(obj) ncol(x = obj@raw_data)))
  obj.names <- lapply(X = obj.list, FUN = function(obj) colnames(x = obj@raw_data))
  if(!all(sapply(X = obj.names[-1], FUN = identical, x = obj.names[[1]]))) stop("Cannot fuse across different cells.", call. = FALSE)
  obj.weights <- do.call(what = rbind, args = lapply(X = obj.list, FUN = function(obj){
    obj.sums <- colSums(x = obj@weight_mat$sample_weight)
  }))
  fused.weight <- matrix(data = 0, nrow = n.obj, ncol = n.cells)
  fused.dist <- list()
  for(i in 1:length(x = obj.list)){
    obj.temp <- obj.list[[i]]
    obj.weight <- colSums(x = obj.temp@weight_mat$sample_weight)
    if(length(x = obj.temp@sketch_id) > 0){
      if(obj.temp@params$pca_dims == 0){
        obj.weight.project <- ProjectSketch(data = t(obj.temp@normalized_data[, obj.temp@sketch_id]), data_query = t(obj.temp@normalized_data[, -obj.temp@sketch_id]), k = project_k, weights_data = obj.weight)
      } else {
        obj.weight.project <- ProjectSketch(data = obj.temp@pca$u[obj.temp@sketch_id, ], data_query = obj.temp@pca$u[-obj.temp@sketch_id, ], k = project_k, weights_data = obj.weight)
      }
      fused.weight[i, obj.temp@sketch_id] <- obj.weight
      fused.weight[i, -obj.temp@sketch_id] <- obj.weight.project
    } else {
      fused.weight[i, ] <- obj.weight
    }
    fused.dist[[i]] <- obj.temp@dist_null
  }
  fused.weight <- apply(X = fused.weight, MARGIN = 2, FUN = function(x) 1-x/sum(x))
  for(i in 1:length(x = obj.list)){
    fused.dist[[i]] <- fused.dist[[i]] * fused.weight[i, ]
    fused.dist[[i]] <- fused.dist[[i]] + t(x = fused.dist[[i]])
  }
  fused.dist <- Reduce(f = "+", x = fused.dist)
  if(nnzero(x = fused.dist)/length(x = fused.dist) < (1 - zero_percent)) fused.dist <- Matrix(data = fused.dist, sparse = TRUE)
  rownames(x = fused.weight) <- unlist(x = lapply(X = obj.list, FUN = function(x) x@project_name))
  colnames(x = fused.weight) <- colnames(x = fused.dist)
  return(list(fused_weight = fused.weight, fused_dist = fused.dist))
}
