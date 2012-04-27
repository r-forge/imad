\name{RADCAL}
\alias{RADCAL}
\title{RADCAL}
\usage{
  RADCAL(inDataSet1, inDataSet2, chisqr_raster,
    noChangeProbThresh = 0.95, minNoChangePixels = NA,
    graph_only = TRUE, return_gains_and_offsets = TRUE,
    apply_to_raster = inDataSet2)
}
\arguments{
  \item{inDataSet1}{A Raster* object of the reference
  image.}

  \item{inDataSet2}{A Raster* object of the image to be
  normalized.}

  \item{chisqr_raster}{Raster. The chi-square image
  generated from iMad (e.g. raster(imad_output,layer=1))}

  \item{noChangeProbThresh}{Numeric. The probability
  threshold (0 <= noChangeProbThresh <= 1) for determining
  no change (default = 0.95).}

  \item{minNoChangePixels}{Logical. NOT SUPPORTED.}

  \item{graph_only}{Logical. NOT SUPPORTED.}

  \item{return_gains_and_offsets}{Logical. Return the gains
  and offsets as a matrix?}

  \item{apply_to_raster}{Raster*. The Raster* to apply the
  gains and the offsets to (default: inDataSet2). Set to NA
  if you don't want to apply the gains and offsets at this
  stage.}
}
\value{
  List (if return_gains_and_offsets==TRUE &&
  !is.na(apply_to_raster)) of the gains and offsets and
  normalized raster. Matrix of gains and rasters if
  !is.na(apply_to_raster). Raster* of the normalized image
  if return_gains_and_offsets==FALSE.
}
\description{
  Perform automated radiometric normalization between two
  rasters following an iMAD calculations.
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

