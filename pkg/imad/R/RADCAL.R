#' RADCAL
#' 
#' Perform automated radiometric normalization between two rasters following an iMAD calculations.
#' 
#' @param inDataSet1 A Raster* object of the reference image.
#' @param inDataSet2 A Raster* object of the image to be normalized.
#' @param chisqr_raster Raster. The chi-square image generated from iMad (e.g. raster(imad_output,layer=1))
#' @param noChangeProbThresh Numeric. The probability threshold (0 <= noChangeProbThresh <= 1) for determining no change (default = 0.95).
#' @param minNoChangePixels Logical. NOT SUPPORTED.
#' @param graph_only Logical. NOT SUPPORTED.
#' @param return_gains_and_offsets Logical. Return the gains and offsets as a matrix?
#' @param apply_to_raster Raster*. The Raster* to apply the gains and the offsets to (default: inDataSet2). Set to NA if you don't want to apply the gains and offsets at this stage.
#' @return List (if return_gains_and_offsets==TRUE && !is.na(apply_to_raster)) of the gains and offsets and normalized raster.  
#' Matrix of gains and rasters if !is.na(apply_to_raster). Raster* of the normalized image if return_gains_and_offsets==FALSE.
#' @author Mort Canty (original code) and Jonathan A. Greenberg (R port).
# @seealso 
# @keywords 
#' @references
#' \itemize{
#' \item {Canty, M.J. and A.A. Nielsen. 2008. Automatic radiometric normalization of multitemporal satellite imagery with the iteratively re-weighted MAD transformation. Remote Sensing of Environment 112:1025-1036.}
#' \item {Nielsen, A.A. 2007. The regularized iteratively reweighted MAD method for change detection in multi- and hyperspectral data. IEEE Transactions on Image Processing 16(2):463-478.}
#' }
# @examples
# 
# \dontrun{
# } 
#' @export

# TODO: Return mask used.
# TODO: Allow for adaptive noChangeProbThresh.
# TODO: Graphing?
# TODO: Return error of gains and offsets.

RADCAL <- function(inDataSet1,inDataSet2,chisqr_raster,noChangeProbThresh=0.95,minNoChangePixels=NA,graph_only=TRUE,
		return_gains_and_offsets=TRUE,apply_to_raster=inDataSet2)
{
	# Do some error checks...
	
	bands=nlayers(inDataSet1)
	
	chisqr_raster_cdf=1-calc(chisqr_raster,fun=function(x) { pchisq(x,bands) })
	pif_raster = chisqr_raster_cdf > noChangeProbThresh
	pif_raster[pif_raster==0] <- NA
	
	pif_N = ncell(pif_raster)- cellStats(pif_raster,'countNA')
	
	# Check for minNoChangePixels here
	
	# Extract all pifs (BE CAREFUL, THIS IS NOT MEMORY SAFE).
	inDataSet1_pifs=inDataSet1[!is.na(pif_raster)]
	inDataSet2_pifs=inDataSet2[!is.na(pif_raster)]
	
	# Some plotting
	# plot(as.matrix(inDataSet1)[,1],as.matrix(inDataSet2)[,1])
	# plot(inDataSet1_pifs[,1]~inDataSet2_pifs[,1])
	
	gains=vector(mode="numeric",length=bands)
	offsets=vector(mode="numeric",length=bands)
	
	# The orthogonal regression may run into sample size issues.
	# Using http://zoonek2.free.fr/UNIX/48_R/09.html
	for(i in 1:bands)
	{	# Orthogonal?
		temp_princomp=princomp(cbind(inDataSet1_pifs[,i],inDataSet2_pifs[,i]))
		gains[i]=temp_princomp$loadings[2,1]/temp_princomp$loadings[1,1]
		offsets[i]=temp_princomp$center[2] - offsets[i] * temp_princomp$center[1]
	}
	
	if(return_gains_and_offsets)
	{
		gains_and_offsets=cbind(gains,offsets)
	}
	
	# Do an error check here:
	return(gains*inDataSet2+offsets)
}