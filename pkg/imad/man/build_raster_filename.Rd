\name{build_raster_filename}
\alias{build_raster_filename}
\title{Build raster filename}
\usage{
  build_raster_filename(filename_base, format = "raster")
}
\arguments{
  \item{filename_base}{Character. The "output filename"
  used by writeRaster (the base filename).}

  \item{format}{Character. The output file type.  Default
  is "raster".}
}
\value{
  Returns a character string that should be usable by
  raster(), brick(), or stack()
}
\description{
  Returns the proper filename to use with raster given the
  base filename and the file format.
}
\author{
  Jonathan A. Greenberg and Robert Hijimans
}
\seealso{
  \code{\link[raster]{writeFormats}},
  \code{\link[raster]{writeRaster}}
}

