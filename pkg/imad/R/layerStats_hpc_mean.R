#' @export

layerStats_hpc_mean <- function(x,na.rm=FALSE, enable_snow=FALSE, cl=NULL, m=2,verbose=FALSE)
{
	require("raster")
	
	if(enable_snow) { 
		require("snowfall")
		if (is.null(cl)) {
			cl <- getCluster()
			on.exit( returnCluster() )
		}
	}
	
	if(verbose) { print("Setting up cluster...")}
	
	if(!enable_snow)
	{
		nodes <- 1
	} else
	{
		if (is.null(cl)) {
			cl <- getCluster()
			on.exit( returnCluster() )
		}
		nodes <- length(cl)
	}
	
	if(verbose) { print("Determining optimal block size...")}
	m <- max(1, round(m))
	tr <- blockSize(x, minblocks=nodes*m )
	if (tr$n < nodes) {
		nodes <- tr$n
	}
	tr$row2 <- tr$row + tr$nrows - 1
	
	i=1:tr$n
	
	nlayers_x=nlayers(x)
	
	if(enable_snow)
	{
		if(verbose) { print("Starting the cluster function...")}
		sumCount <- clusterMap(cl,function(i,x,tr,na.rm) 
				{
					r <- getValues(crop(x, extent(x, r1=tr$row[i], r2=tr$row2[i], c1=1, c2=ncol(x))))
					if(class(r)=="matrix")
					{
						temp_sum=colSums(as.array(r),na.rm=na.rm)
						temp_count=colSums(!is.na(as.array(r)))
					} else
					{
						temp_sum=sum(r,na.rm=na.rm)
						temp_count=colSums(!is.na(as.array(r)))
					}
					return(cbind(temp_sum,temp_count))
				},
				i,MoreArgs=list(x=x,tr=tr,na.rm=na.rm))
		if(nlayers_x > 1){
			layersums=rowSums(sapply(sumCount,function(x) { x[,1] }),na.rm=na.rm)
			layercounts=rowSums(sapply(sumCount,function(x) { x[,2] }),na.rm=na.rm)
		} else
		{
			layersums=sum(sapply(sumCount,function(x) { x[,1] }),na.rm=na.rm)
			layercounts=sum(sapply(sumCount,function(x) { x[,2] }),na.rm=na.rm)
		}
		return(layersums/layercounts)
	} else
	{
		# Put in a safe version of sum here
		return(NULL)
	}
}