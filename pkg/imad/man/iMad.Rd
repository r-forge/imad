\name{iMad}
\alias{iMad}
\title{iMad (optimized port)}
\usage{
  iMad(inDataSet1, inDataSet2, maxiter = 100, lam = 0,
    output_basename, verbose = FALSE,
    auto_extract_overlap = TRUE,
    reuse_existing_overlap = TRUE, corr_thresh = 0.001,
    delta = 0.001, format = "raster", mask1, mask2,
    force_extent = TRUE, ...)
}
\arguments{
  \item{inDataSet1}{A Raster* object of the first image.}

  \item{inDataSet2}{A Raster* object of the second image.}

  \item{maxiter}{Numeric (>= 1). The maximum number of
  iterations.}

  \item{lam}{Numeric. The penalization function.  CURRENTLY
  UNSUPPORTED.}

  \item{output_basename}{Character. The basename (including
  path) for the output files.}

  \item{verbose}{Logical. Print out debugging information?}

  \item{auto_extract_overlap}{Logical. Extract the overlap
  zones between the images?}

  \item{reuse_existing_overlap}{Logical. If the algorithm
  detects pre-create overlaps, use them?}

  \item{delta}{Numeric. The smallest change in canonical
  correlates to end the program.}

  \item{corr_thresh}{Numeric. Used for situations where the
  canonical correlates are all nearly 1.0 (how close to 1.0
  does it need to be to stop).}

  \item{format}{Character. The output format of the rasters
  (see ?writeFormats).}

  \item{mask1}{A Raster* object representing a mask to be
  used for inDataSet1.}

  \item{mask2}{A Raster* object representing a mask to be
  used for inDataSet2.}

  \item{force_extent}{Logical. Attempt to force the input
  files (and masks) to be the same extent?}

  \item{...}{Passed to various raster functions (see
  writeRaster). Important ones include format= and
  overwrite=TRUE/FALSE.  datatype should be left as 'FLT4S'
  for proper functioning.}
}
\value{
  Returns a RasterStack object where the first layer is the
  chisquare image, and the subsequent layers are the iMad
  layers.
}
\description{
  Perform iteratively re-weighted multivariate alteration
  detection.
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
  \code{\link[raster]{writeRaster}}
}

