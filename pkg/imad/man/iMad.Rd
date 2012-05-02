\name{iMad}
\alias{iMad}
\title{iMad (optimized port)}
\usage{
  iMad(inDataSet1, inDataSet2, pos, mask1, mask2,
    mask1_band = 1, mask2_band = 1, output_basename,
    format = "raster", maxiter = 100, lam = 0,
    delta = 0.001, corr_thresh = 0.001,
    auto_extract_overlap = TRUE,
    reuse_existing_raster = TRUE, force_extent = TRUE,
    use_clusterR = FALSE, cl = NULL, verbose = FALSE, ...)
}
\arguments{
  \item{inDataSet1}{A Raster* object of the first image.}

  \item{inDataSet2}{A Raster* object of the second image.}

  \item{pos}{Integer vector.  A vector of bands to use from
  each image.  Default is use all bands.}

  \item{mask1}{(Optional) A Raster* object representing a
  mask to be used for inDataSet1 or a numeric value to be
  used as the mask value.}

  \item{mask2}{(Optional) A Raster* object representing a
  mask to be used for inDataSet2 or a numeric value to be
  used as the mask value.}

  \item{mask1_band}{(Optional) The band from inDataSet1 to
  use for masking (only if class(mask1)=="numeric").}

  \item{mask2_band}{(Optional) The band from inDataSet2 to
  use for masking (only if class(mask1)=="numeric").}

  \item{output_basename}{Character. The basename (including
  path) for the output files.}

  \item{format}{Character. The output format of the rasters
  (see ?writeFormats).  Default is "raster".}

  \item{maxiter}{Numeric (>= 1). The maximum number of
  iterations.  Default is 100.}

  \item{lam}{Numeric. The penalization function.  CURRENTLY
  UNSUPPORTED.}

  \item{corr_thresh}{Numeric. Used for situations where the
  canonical correlates are all nearly 1.0 (how close to 1.0
  does it need to be to stop).}

  \item{delta}{Numeric. The smallest change in canonical
  correlates to end the program.}

  \item{auto_extract_overlap}{Logical. Extract the overlap
  zones between the images?}

  \item{reuse_existing_raster}{Logical. If the algorithm
  detects pre-create overlaps, use them?}

  \item{force_extent}{Logical. Attempt to force the input
  files (and masks) to be the same extent?}

  \item{use_clusterR}{Logical. Use clusterR to
  (potentially) speed up calculations on a cluster?
  Default=FALSE.}

  \item{cl}{Cluster.  If not assigned, the program will
  attempt to figure it out.}

  \item{verbose}{Logical. Print out debugging information?}

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

