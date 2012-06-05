# Jonathan Greenberg and Robert Hijmans
# Date : April 2012
# Version 1.0
# Licence GPL v3

# Computation of the weighted covariance and (optionally) weighted means of bands in an Raster.
# based on code by Mort Canty

#' @export

layerStats_hpc <- function(x, stat, w, asSample=FALSE, na.rm=FALSE, enable_snow=TRUE, cl=NULL, m=2, verbose=FALSE, ...) {
	if(enable_snow) { 
		require("snowfall")
		if (is.null(cl)) {
			cl <- getCluster()
			on.exit( returnCluster() )
		}
	}
	
	stat <- tolower(stat)
	stopifnot(stat %in% c('weighted.cov', 'sum','mean'))
	stopifnot(is.logical(asSample) & !is.na(asSample))
	
	if (stat == 'weighted.cov') {
		return(layerStats_hpc_weighted.cov(x,w,na.rm=na.rm, enable_snow=enable_snow, cl=cl, m=m,verbose=verbose))
	}
	
	if(stat=="sum")
	{
		return(layerStats_hpc_sum(x,na.rm=na.rm, enable_snow=enable_snow, cl=cl, m=m,verbose=verbose))
	}
	
	if(stat=="mean")
	{
		return(layerStats_hpc_mean(x,na.rm=na.rm, enable_snow=enable_snow, cl=cl, m=m,verbose=verbose))
	}
	
	
	
}


