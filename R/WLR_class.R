library(Matrix)
setClassUnion(name = 'matrice', members = c("matrix", "dgCMatrix", "dsCMatrix"))

#' The WLR S4 Class
#'
#' The WLR object with the slot information listed as follow:
#'
#' @slot project_name Name of the project.
#' @slot raw_data Raw data. A d x M matrix with d rows of features and M columns of data points.
#' @slot normalized_data Normalized data. Same shape as raw_data.
#' @slot scaled_data Scaled data (z-score). Same shape as raw_data.
#' @slot pca Principal component analysis result, see \code{\link[irlba]{irlba}}.
#' @slot dist_null Null nearest neighbor M x M distance matrix.
#' @slot sketch_id Geomertric sketching cell IDs.
#' @slot sketch_dist Geomertric sketching distance matrix.
#' @slot weight_mat List of feature and sample weight matrice:
#' \itemize{
#'   \item feature_weight, n x d binary matrix with n rows of bootstrap iterations and d columns of features where 0 means feature not sampled and 1 means sampled.
#'   \item sample_weight, n x M matrix with n rows of bootstrap iterations and M columns of data points where each entry represents weight.
#'   \item perturb_mat, d x M matrix with d rows of features and M columns of data points where each entry represents the relative importance of a feature to a data point.
#'   }
#' @slot dist_mat Permuted distance matrix.
#' @slot params Commands used.
#' @name WLR-class
#' @rdname WLR-class
#' @exportClass WLR
#'
setClass(Class = "WLR",
         slots = c(
           project_name = "character",
           raw_data = "matrice",
           normalized_data = "matrice",
           scaled_data = "matrix",
           pca = "list",
           dist_null = "matrice",
           sketch_id = "numeric",
           sketch_dist = "matrice",
           weight_mat = "list",
           dist_mat = "matrice",
           params = "list"
         )
)

setMethod(f = "show",signature = "WLR",
          definition = function(object){
            cat(object@project_name, "WLR object", "\n")
            cat(nrow(object@raw_data), "features across", ncol(object@raw_data), "samples", "\n")
          }
)
