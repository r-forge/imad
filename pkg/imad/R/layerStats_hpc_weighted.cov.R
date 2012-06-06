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
			if(verbose) { print("Calculating weight sum...")}
			sumw=layerStats_hpc(w,'sum',na.rm=na.rm)
			if(verbose) { print("Multiplying x by w...")}
			xw = calc_hpc(stack(w,x),
					fun=function(x) 
					{ 
						nlayers_x=nlayers(x)
						w=raster(x,layer=1)
						pos=2:nlayers_x
						x_image=stack(mapply(function(band,inbrick){raster(inbrick,layer=band)},band=pos,MoreArgs=list(inbrick=x)))
						return(w*x_image)
					})	
			if(verbose) { print("Calculated weighted means...")}
			means=layerStats_hpc(xw,'sum',na.rm=na.rm)/sumw
#			sumw <- sumw - asSample
			if(verbose) { print("sqrt(wt)...")}
#			w_sqrt = calc_hpc(x=w,fun=sqrt)
#			x = stack(clusterMap(cl,fun=function(x,means,w_sqrt) { (x - means) * w_sqrt },
#							x=raster_to_list(x),MoreArgs=list(means=means,w_sqrt=w_sqrt)))
			if(verbose) { print("(x-means)*w_sqrt...")}
			x=calc_hpc(stack(w,x),args=list(means=means),
				fun=function(x,means)
				{
					nlayers_x=nlayers(x)
					pos=2:nlayers_x
					w=raster(x,layer=1)
					x_image=stack(mapply(function(band,inbrick){raster(inbrick,layer=band)},band=pos,MoreArgs=list(inbrick=x)))
					out=(x_image-means)*calc(w,sqrt)
					print(out)
					return(out)
				},verbose=verbose
			)
			
		} else
		{
			sumw <- cellStats(w, stat='sum', na.rm=na.rm) 
			xw=x*w
			means <- cellStats(xw, stat='sum', na.rm=na.rm) / sumw
			sumw <- sumw - asSample
			w_sqrt=sqrt(w)
			x <- (x - means) * w_sqrt
		}
		
		if(verbose) { print("creating ij") }
		ij=which(upper.tri(matrix(ncol=nl,nrow=nl),diag=TRUE),arr.ind=TRUE)
		ij_idx=as.list(1:(dim(ij)[1]))
		ij_list=mapply(function(ij_idx,ij) { ij[ij_idx,] },ij_idx=ij_idx,MoreArgs=list(ij=ij),SIMPLIFY=FALSE)
		
		if(enable_snow) { 
		
#			v_list=clusterMap(cl,fun=function(ij,x,na.rm,sumw) { 
#						i <- ij[1]
#						j <- ij[2]
#						r <- raster(x,layer=i)*raster(x,layer=j)
#						v <- cellStats(r, stat='sum', na.rm=na.rm) / sumw
#						return(v)
#					},
#					ij=ij_list,MoreArgs=list(x=x,na.rm=na.rm,sumw=sumw))
	# Need to spawn mini clusters for this, or just let it go sequentially...
			if(verbose) { print("creating v_list") }
			v_list=mapply(FUN=function(ij,x,na.rm,sumw) { 
						i <- ij[1]
						j <- ij[2]
						rasteri=raster(x,layer=i)
						rasterj=raster(x,layer=j)
#						r <- raster(x,layer=i)*raster(x,layer=j)
						r=calc_hpc(stack(rasteri,rasterj),
								fun=function(x)
								{
									out=raster(x,layer=1)*raster(x,layer=2)
									return(out)
								})
#						v <- cellStats(r, stat='sum', na.rm=na.rm) / sumw
						v=layerStats_hpc(r,stat='sum',na.rm=na.rm)/sumw
						return(v)
					},
					ij=ij_list,MoreArgs=list(x=x,na.rm=na.rm,sumw=sumw))	
			if(verbose) { print("done with v_list") }
			for(k in 1:(length(v_list)))
			{
				i=ij_list[[k]][1]
				j=ij_list[[k]][2]
				mat[j,i] <- mat[i,j] <- v_list[[k]]
			}
			if(verbose) { print("Done making the matrix...")}
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