\name{cov.wt.raster}
\alias{cov.wt.raster}
\title{Weighted Covariance Matrices and Means for Rasters}
\usage{
  cov.wt.raster(dm, wt, return_means = TRUE)
}
\arguments{
  \item{dm}{A multiband Raster* to use in calculating
  weighted covariance matrices.}

  \item{wt}{A Raster* object of the weights (should have
  the same extent as dm).}

  \item{return_means}{Logical. Return the weighted per-band
  means?}
}
\value{
  Returns either a covariance matrix of dimensions
  nlayers(dm) x nlayers(dm) or a list containing the
  covariance matrix and the weighted per-band means.
}
\description{
  Calculates the weighted covariance and (optionally)
  weighted means of bands in an Raster.
}
\author{
  Mort Canty (original code) and Jonathan A. Greenberg (R
  port).
}
\references{
  \itemize{ \item {Canty, M.J. and A.A. Nielsen. 2008.
  Automatic radiometric normalization of multitemporal
  satellite imagery with the iteratively re-weighted MAD
  transformation. Remote Sensing of Environment
  112:1025-1036.} \item {Nielsen, A.A. 2007. The
  regularized iteratively reweighted MAD method for change
  detection in multi- and hyperspectral data. IEEE
  Transactions on Image Processing 16(2):463-478.} }
}
\seealso{
  \code{\link{cov.wt}}, \code{\link{weighted.mean}}
}

