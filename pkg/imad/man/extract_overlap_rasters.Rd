\name{extract_overlap_rasters}
\alias{extract_overlap_rasters}
\title{extract_overlap_rasters}
\usage{
  extract_overlap_rasters(raster1, raster2, filename1,
    filename2, raster1_crop = TRUE, raster2_crop = TRUE,
    verbose = TRUE, ...)
}
\arguments{
  \item{raster1}{The first raster to extract the overlap
  zone from.}

  \item{raster2}{The second raster to extract the overlap
  zone from.}

  \item{filename1}{The output filename of raster1
  (optional).}

  \item{filename2}{The output filename of raster2
  (optional).}

  \item{raster1_crop}{Logical. Proceed with cropping and
  returning raster1?}

  \item{raster2_crop}{Logical. Proceed with cropping and
  returning raster2?}

  \item{verbose}{Logical. Provide debugging feedback?}

  \item{...}{Other parameters to pass to crop.}
}
\description{
  Coerces two rasters to have a shared subset of their
  common overlap.
}
\author{
  Jonathan A. Greenberg \email{STARStools@estarcion.net}
}
\seealso{
  \code{\link[raster]{crop}}
}

