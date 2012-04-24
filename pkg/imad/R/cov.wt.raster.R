#' Weighted Covariance Matrices and Means for Rasters
#' 
#' Calculates the weighted covariance and (optionally) weighted means of bands in an Raster.
#' 
#' @param dm A multiband Raster* to use in calculating weighted covariance matrices.
#' @param wt A Raster* object of the weights (should have the same extent as dm).
#' @param return_means Logical. Return the weighted per-band means?
#' @return Returns either a covariance matrix of dimensions nlayers(dm) x nlayers(dm) or a list containing the covariance matrix and the weighted per-band means.
#' @author Mort Canty (original code) and Jonathan A. Greenberg (R port).
#' @seealso \code{\link{cov.wt}}, \code{\link{weighted.mean}}
#' @references
#' \itemize{
#' \item {Canty, M.J. and A.A. Nielsen. 2008. Automatic radiometric normalization of multitemporal satellite imagery with the iteratively re-weighted MAD transformation. Remote Sensing of Environment 112:1025-1036.}
#' \item {Nielsen, A.A. 2007. The regularized iteratively reweighted MAD method for change detection in multi- and hyperspectral data. IEEE Transactions on Image Processing 16(2):463-478.}
#' }
# @keywords {weighted covariance matrix}
# {weighted means}
# @examples
# \dontrun{
# } 
#' @export

cov.wt.raster <- function(dm,wt,return_means=TRUE)
{
	if(missing(wt))
	{
		wt=raster(dm,layer=1)*0+1
	}
	
	N = nlayers(dm)
	n = ncell(wt)
	sumw = cellStats(wt,stat='sum')
	layer=as.list(1:N)
	means = mapply(function(wt,dm,sumw,layer) { cellStats(wt*raster(dm,layer=layer)/sumw,stat='sum',na.rm=TRUE) },
			layer,MoreArgs=list(wt=wt,dm=dm,sumw=sumw))
	dmc = (dm - means)*sqrt(wt)
	
	# We should do this more efficiently but...
	covmat=matrix(nrow=N,ncol=N)
	for(i in 1:N)
	{
		for(j in 1:N)
		{
			covmat[i,j]=cellStats(raster(dmc,layer=i)*raster(dmc,layer=j),stat='sum',na.rm=TRUE)/sumw	
		}
	}
	if(return_means)
	{
		cov.wt.raster <- list(covmat,means)
		names(cov.wt.raster) <- c("covariance","means")
		return(cov.wt.raster)
	} else
	{
		return(covmat)
	}
}