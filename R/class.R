#' @importClassesFrom Matrix dgCMatrix dsCMatrix
#' 
NULL

setClassUnion(name = 'matrices', members = c("matrix", "dgCMatrix", "dsCMatrix"))

#' The FuseNet S4 Class
#'
#' The FuseNet object with the slot information listed as follow:
#'
#' @slot project_name Name of the project.
#' @slot raw_data Raw data. A d x M matrix with d rows of features and M columns of data points.
#' @slot normalized_data Normalized data. Same shape as raw_data.
#' @slot scaled_data Scaled data (z-score). Same shape as raw_data.
#' @slot pca Principal component analysis result, see \code{\link[irlba]{irlba}}.
#' @slot dist_null Null nearest neighbor M x M distance matrix.
#' @slot sketch_id Geomertric sketching cell IDs.
#' @slot sketch_dist Geomertric sketching distance matrix.
#' @slot weight_mat List of feature and sample weight matrices:
#' \itemize{
#'   \item feature_weight, n x d binary matrix with n rows of bootstrap iterations and d columns of features where 0 means feature not sampled and 1 means sampled.
#'   \item sample_weight, n x M matrix with n rows of bootstrap iterations and M columns of data points where each entry represents weight.
#'   \item perturb_mat, d x M matrix with d rows of features and M columns of data points where each entry represents the relative importance of a feature to a data point.
#'   }
#' @slot dist_mat Permuted distance matrix.
#' @slot params Commands used.
#' @name FuseNet-class
#' @rdname FuseNet-class
#' @exportClass FuseNet
#'
setClass(Class = "FuseNet",
         slots = c(
           project_name = "character",
           raw_data = "matrices",
           normalized_data = "matrices",
           scaled_data = "matrix",
           pca = "list",
           dist_null = "matrices",
           sketch_id = "numeric",
           sketch_dist = "matrices",
           weight_mat = "list",
           dist_mat = "matrices",
           params = "list"
         )
)

setMethod(f = "show",signature = "FuseNet",
          definition = function(object){
            cat(object@project_name, "FuseNet object", "\n")
            cat(nrow(object@raw_data), "features across", ncol(object@raw_data), "samples", "\n")
          }
)
