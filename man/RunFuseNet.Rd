% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fuse.R
\name{RunFuseNet}
\alias{RunFuseNet}
\title{Run Data Fusion}
\usage{
RunFuseNet(
  object,
  n_iters = 100,
  ratio = 0.05,
  pca_dims = 0,
  k = 100,
  t = 0,
  norm_type = c("l1", "l2"),
  return_perturb_mat = FALSE,
  n_cores = NULL,
  ...
)
}
\arguments{
\item{object}{A FuseNet object.}

\item{n_iters}{Number of bootstrapping iterations. Default is 100.}

\item{ratio}{Fraction of features to be downsampled in the original data matrix. Default is 0.05 aka 5\%.}

\item{pca_dims}{Number of principle components. Default is 0 and PCA is not run.}

\item{k}{Number of nearest neighbors used. Default is 100.}

\item{t}{Matrix power used for the distance matrix. Default is 0 and powering is not performed.}

\item{norm_type}{Type of norm used:
\itemize{
  \item l1, L1-like norm. See details \code{\link[FuseNet]{L1Norm}}.
  \item l2, L1-like norm. See details \code{\link[FuseNet]{L2Norm}}.
  }}

\item{return_perturb_mat}{Whether to return the perturb matrix. Default is FALSE.}

\item{n_cores}{Number of cores used. Default is to use all existing cores. See details \code{\link[parallel]{makeCluster}}.}

\item{...}{Additional parameters pass to \code{\link[parallel]{makeCluster}}.}
}
\value{
Returns a FuseNet object.
}
\description{
Run Data Fusion.
}
\examples{
{
object <- InitiateFuseNet(t(iris[,1:4]), project_name = "FuseNet", k = 3)
object <- RunFuseNet(object, n_iters = 1, k = 10, ratio = 0.5, n_cores = 1)
}

}
\concept{fuse}
