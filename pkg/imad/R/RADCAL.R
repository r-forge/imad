

# TODO: Return mask used.
# TODO: Allow for adaptive noChangeProbThresh based on sample size and 

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