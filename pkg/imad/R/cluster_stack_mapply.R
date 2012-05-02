
cluster_stack_mapply <- function(x, fun, args=NULL, filename='', cl=NULL, m=2, ...)
{
	
	raster_to_list=function(inraster)
	{
		inraster_nlayers=nlayers(inraster)
		pos=1:inraster_nlayers
		return(
			mapply(function(band,inraster){ raster(inraster,layer=band) },
			band=pos,MoreArgs=list(inraster=inraster),SIMPLIFY=FALSE)
		)
		
	}

	if (is.null(cl)) {
		cl <- getCluster()
		on.exit( returnCluster() )
	}
	
	x_list=raster_to_list(x)
	
#	if (!is.null(args)) {
#		stopifnot(is.list(args))
#		
#		clusfun <- function(fun, i) {
##			r <- crop(x, extent(out, r1=tr$row[i], r2=tr$row2[i], c1=1, c2=ncol(out)))
#			r <- do.call(fun, c(r, args))
##			getValues(r)
#		}
#		
#	} else {
#		
#		clusfun <- function(fun, i) {
##			r <- crop(x, extent(out, r1=tr$row[i], r2=tr$row2[i], c1=1, c2=ncol(out)))
#			r <- fun(r)
##			getValues(r)
#		}
#	}
	
	x_out = stack(clusterMap(cl,fun=clusfun,x=x_list,MoreArgs=args))
	return(x_out)
}
