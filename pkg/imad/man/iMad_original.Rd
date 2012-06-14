\name{iMad_original}
\alias{iMad_original}
\title{iMad Original (direct python port)}
\usage{
  iMad_original(inDataSet1, inDataSet2, maxiter = 100,
    lam = 0, output_basename, verbose = FALSE, ...)
}
\arguments{
  \item{inDataSet1}{A Raster* object of the first image.}

  \item{inDataSet2}{A Raster* object of the second image.}

  \item{maxiter}{Numeric (>= 1). The maximum number of
  iterations.}

  \item{lam}{Numeric. The penalization function.}

  \item{output_basename}{Character. The basename (including
  path) for the output files.}
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
\examples{
\dontrun{
}
}
\author{
  Mort Canty (original code) and Jonathan A. Greenberg (R
  port).
}

