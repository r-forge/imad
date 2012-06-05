# Jonathan Greenberg and Robert Hijmans
# Date : April 2012
# Version 1.0
# Licence GPL v3

# Computation of the weighted covariance and (optionally) weighted means of bands in an Raster.
# based on code by Mort Canty

#' @export

layerStats_hpc <- function(x, stat, w, asSample=TRUE, na.rm=FALSE, enable_snow=TRUE, cl=NULL, verbose=FALSE, ...) {
	if(enable_snow) { 
		require("snowfall")
		if (is.null(cl)) {
			cl <- getCluster()
			on.exit( returnCluster() )
		}
	}
	
	stat <- tolower(stat)
	stopifnot(stat %in% c('cov', 'weighted.cov', 'pearson', 'sum','mean'))
	stopifnot(is.logical(asSample) & !is.na(asSample))

	nl <- nlayers(x)
	n <- ncell(x)
	mat <- matrix(NA, nrow=nl, ncol=nl)
	colnames(mat) <- rownames(mat) <- layerNames(x)
	
	if (stat == 'weighted.cov') {
		if (missing(w))	{
			stop('to compute weighted covariance a weights layer should be provided')
		}
		stopifnot( nlayers(w) == 1 )

#		if (na.rm) {
#		# a cell is set to NA if it is NA in any layer. That is not ideal, but easier and quicker
#			nas <- calc(x, function(i) sum(i)) * w
#			x <- mask(x, nas)
#			w <- mask(w, nas)
#		}

		sumw <- cellStats(w, stat='sum', na.rm=na.rm) 
		if(enable_snow)
		{
			xw = calc_hpc(stack(w,x),
				fun=function(x) 
				{ 
					nlayers_x=nlayers(x)
					w=raster(x,layer=1)
					pos=2:nlayers_x
					x_image=stack(mapply(function(band,inbrick){raster(inbrick,layer=band)},band=pos,MoreArgs=list(inbrick=x)))
					return(w*x_image)
				})	
			means <- cellStats(xw, stat='sum', na.rm=na.rm) / sumw
		} else
		{
			xw=x*w
			means <- cellStats(xw, stat='sum', na.rm=na.rm) / sumw
		}
		
		sumw <- sumw - asSample
		
		if(enable_snow) { 
			w_sqrt = calc_hpc(x=w,fun=sqrt)
			x = stack(clusterMap(cl,fun=function(x,means,w_sqrt) { (x - means) * w_sqrt },
				x=raster_to_list(x),MoreArgs=list(means=means,w_sqrt=w_sqrt)))
		} else
		{
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
		
	} else if (stat == 'cov') {

		means <- cellStats(x, stat='mean', na.rm=na.rm) 
		x <- (x - means)
		
		for(i in 1:nl) {
			for(j in i:nl) {
				r <- raster(x, layer=i) * raster(x, layer=j)
				if (na.rm) {
					v <- cellStats(r, stat='sum', na.rm=na.rm) / (n - cellStats(r, stat='countNA') - asSample)
				} else {
					v <- cellStats(r, stat='sum', na.rm=na.rm) / (n - asSample)
				}
				mat[j,i] <- mat[i,j] <- v
			}
		}
		covar <- list(mat, means)
		names(covar) <- c("covariance", "mean")
		return(covar)		
		
	} else if (stat == 'pearson') {

		means <- cellStats(x, stat='mean', na.rm=na.rm) 
		sds <- cellStats(x, stat='sd', na.rm=na.rm) 
		x <- (x - means)
		
		for(i in 1:nl) {
			for(j in i:nl) {
				r <- raster(x, layer=i) * 	raster(x, layer=j)
				if (na.rm) {
					v <- cellStats(r, stat='sum', na.rm=na.rm) / ((n - cellStats(r, stat='countNA') - asSample) * sds[i] * sds[j])
				} else {
					v <- cellStats(r, stat='sum', na.rm=na.rm) / ((n - asSample) * sds[i] * sds[j])
				}
				mat[j,i] <- mat[i,j] <- v
			}
		}
		covar <- list(mat, means)
		names(covar) <- c("pearson correlation coefficient", "mean")
		return(covar)
		
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


