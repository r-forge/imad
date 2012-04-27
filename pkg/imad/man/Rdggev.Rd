\name{Rdggev}
\alias{Rdggev}
\title{Rdggev: R wrapper for LAPACK's dggev}
\usage{
  Rdggev(A, B, JOBVL = FALSE, JOBVR = TRUE)
}
\arguments{
  \item{A}{NxN nonsymmetric matrix.}

  \item{B}{NxN nonsymmetric matrix.}

  \item{JOBVL}{Logical. Compute the left generalized
  eigenvectors?}

  \item{JOBVR}{Logical. Compute the right generalized
  eigenvectors?}
}
\value{
  List of input/output parameters of dggev as well as a
  calculation of the eigenvalues (GENEIGENVALUES) (see
  caveats).
}
\description{
  Computes for a pair of N-by-N real nonsymmetric matrices
  (A,B) the generalized eigenvalues, and optionally, the
  left and/or right generalized eigenvectors.
}
\note{
  Please see this thread for a discussion about this
  function:
  \url{http://tolstoy.newcastle.edu.au/R/e17/help/12/04/10951.html}
}
\examples{
A <- matrix(c(1457.738, 1053.181, 1256.953,
             1053.181, 1213.728, 1302.838,
             1256.953, 1302.838, 1428.269), nrow=3, byrow=TRUE)
B <- matrix(c(4806.033, 1767.480, 2622.744,
             1767.480, 3353.603, 3259.680,
             2622.744, 3259.680, 3476.790), nrow=3, byrow=TRUE)
rdggev_out <- Rdggev(A,B,JOBVL=FALSE,JOBVR=TRUE)
print("Generalized eigenvalues:")
rdggev_out$GENEIGENVALUES
print("Right eigenvectors:")
rdggev_out$VR
}
\author{
  Jonathan A. Greenberg
}
\references{
  \url{http://www.netlib.org/lapack/double/dggev.f}
}

