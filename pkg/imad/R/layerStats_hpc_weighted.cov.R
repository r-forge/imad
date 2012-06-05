#' @export

layerStats_hpc_weighted.cov <- function(x,w,na.rm=FALSE, asSample=FALSE,enable_snow=FALSE, cl=NULL, m=2,verbose=FALSE)
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
	
	nl <- nlayers(x)
	n <- ncell(x)
	mat <- matrix(NA, nrow=nl, ncol=nl)
	colnames(mat) <- rownames(mat) <- names(x)
	
	i=1:tr$n
	
		if (missing(w))	{
			stop('to compute weighted covariance a weights layer should be provided')
		}
		stopifnot( nlayers(w) == 1 )
		
		if(enable_snow)
		{
			sumw=layerStats_hpc(w,'sum',na.rm=na.rm)
			xw = calc_hpc(stack(w,x),
					fun=function(x) 
					{ 
						nlayers_x=nlayers(x)
						w=raster(x,layer=1)
						pos=2:nlayers_x
						x_image=stack(mapply(function(band,inbrick){raster(inbrick,layer=band)},band=pos,MoreArgs=list(inbrick=x)))
						return(w*x_image)
					})	
			means=layerStats_hpc(xw,'sum',na.rm=na.rm)/sumw
#			sumw <- sumw - asSample
			w_sqrt = calc_hpc(x=w,fun=sqrt)
			x = stack(clusterMap(cl,fun=function(x,means,w_sqrt) { (x - means) * w_sqrt },
							x=raster_to_list(x),MoreArgs=list(means=means,w_sqrt=w_sqrt)))
		} else
		{
			sumw <- cellStats(w, stat='sum', na.rm=na.rm) 
			xw=x*w
			means <- cellStats(xw, stat='sum', na.rm=na.rm) / sumw
			sumw <- sumw - asSample
			w_sqrt=sqrt(w)
			x <- (x - means) * w_sqrt
		}
		
		ij=which(upper.tri(matrix(ncol=nl,nrow=nl),diag=TRUE),arr.ind=TRUE)
		ij_idx=as.list(1:(dim(ij)[1]))
		ij_list=mapply(function(ij_idx,ij) { ij[ij_idx,] },ij_idx=ij_idx,MoreArgs=list(ij=ij),SIMPLIFY=FALSE)
		
		if(enable_snow) { 
			v_list=clusterMap(cl,fun=function(ij,x,na.rm,sumw) { 
						i <- ij[1]
						j <- ij[2]
						r <- raster(x,layer=i)*raster(x,layer=j)
						v <- cellStats(r, stat='sum', na.rm=na.rm) / sumw
						return(v)
					},
					ij=ij_list,MoreArgs=list(x=x,na.rm=na.rm,sumw=sumw))
			for(k in 1:(length(v_list)))
			{
				i=ij_list[[k]][1]
				j=ij_list[[k]][2]
				mat[j,i] <- mat[i,j] <- v_list[[k]]
			}
		} else
		{
			for(i in 1:nl) {
				for(j in i:nl) {
					r <- raster(x, layer=i) * raster(x,layer=j)
					v <- cellStats(r, stat='sum', na.rm=na.rm) / sumw
					mat[j,i] <- mat[i,j] <- v
					
				}
			}
		}
		cov.w <- list(mat, means)
		names(cov.w) <- c("weighted covariance", "weighted mean")
		return(cov.w)		
}