#' Rdggev: R wrapper for LAPACK's dggev
#' 
#' Computes for a pair of N-by-N real nonsymmetric matrices (A,B) the generalized eigenvalues, and optionally, the left and/or right generalized eigenvectors.
#' 
#' @param A NxN nonsymmetric matrix.
#' @param B NxN nonsymmetric matrix.
#' @param JOBVL Logical. Compute the left generalized eigenvectors?
#' @param JOBVR Logical. Compute the right generalized eigenvectors?
#' @return List of input/output parameters of dggev as well as a calculation of the eigenvalues (GENEIGENVALUES) (see caveats).
#' @author Jonathan A. Greenberg
# @seealso 
# @keywords 
#' @note Please see this thread for a discussion about this function: \url{http://tolstoy.newcastle.edu.au/R/e17/help/12/04/10951.html}
#' @references \url{http://www.netlib.org/lapack/double/dggev.f}
#' @examples
#' A <- matrix(c(1457.738, 1053.181, 1256.953,
#'              1053.181, 1213.728, 1302.838,
#'              1256.953, 1302.838, 1428.269), nrow=3, byrow=TRUE)
#' B <- matrix(c(4806.033, 1767.480, 2622.744,
#'              1767.480, 3353.603, 3259.680,
#'              2622.744, 3259.680, 3476.790), nrow=3, byrow=TRUE)
#' rdggev_out <- Rdggev(A,B,JOBVL=FALSE,JOBVR=TRUE)
#' print("Generalized eigenvalues:")
#' rdggev_out$GENEIGENVALUES
#' print("Right eigenvectors:")
#' rdggev_out$VR
#' @export

Rdggev <- function(A,B,JOBVL=FALSE,JOBVR=TRUE)
{
	# R implementation of the DGGEV LAPACK function (with generalized eigenvalue computation)
	# See http://www.netlib.org/lapack/double/dggev.f
	
	# coded by Jonathan A. Greenberg <imad@estarcion.net>
	# Contributions from Berend Hasselman.
	
	if( .Platform$OS.type == "windows" ) {
		Lapack.so <- file.path(R.home("bin"),paste("Rlapack",.Platform$dynlib.ext,sep=""))
	} else {
		Lapack.so <- file.path(R.home("modules"),paste("lapack",.Platform$dynlib.ext,sep=""))
	}
	
	dyn.load(Lapack.so)
	
	if(JOBVL)
	{
		JOBVL="V"
	} else
	{
		JOBVL="N"
	}
	
	if(JOBVR)
	{
		JOBVR="V"
	} else
	{
		JOBVR="N"
	}
	
	if(!is.matrix(A)) stop("Argument A should be a matrix")
	if(!is.matrix(B)) stop("Argument B should be a matrix")
	dimA <- dim(A)
	if(dimA[1]!=dimA[2]) stop("A must be a square matrix")
	dimB <- dim(B)
	if(dimB[1]!=dimB[2]) stop("B must be a square matrix")
	if(dimA[1]!=dimB[1]) stop("A and B must have the same dimensions")
	
	if( is.complex(A) ) stop("A may not be complex")
	if( is.complex(B) ) stop("B may not be complex")
	
	# Input parameters
	N=dim(A)[[1]]
	LDA=N
	LDB=N
	LDVL=N
	LDVR=N
	LWORK=as.integer(max(1,8*N))
	
	Rdggev_out <- .Fortran("dggev", JOBVL, JOBVR, N, A, LDA, B, LDB, double(N), double(N), double(N),
			array(data=0,dim=c(LDVL,N)), LDVL, array(data=0,dim=c(LDVR,N)), LDVR, double(max(1,LWORK)), LWORK, integer(1))
	
	names(Rdggev_out)=c("JOBVL","JOBVR","N","A","LDA","B","LDB","ALPHAR","ALPHAI",
			"BETA","VL","LDVL","VR","LDVR","WORK","LWORK","INFO")
	
	# simplistic calculation of eigenvalues (see caveat in http://www.netlib.org/lapack/double/dggev.f)
	if( all(Rdggev_out$ALPHAI==0) )
		Rdggev_out$GENEIGENVALUES <- Rdggev_out$ALPHAR/Rdggev_out$BETA
	else
		Rdggev_out$GENEIGENVALUES <- complex(real=Rdggev_out$ALPHAR, imaginary=Rdggev_out$ALPHAI)/Rdggev_out$BETA  
	
	return(Rdggev_out)
}