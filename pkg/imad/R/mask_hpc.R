#' @export

mask_hpc <- function(x, mask, cl=NULL, disable_cl=FALSE,inmemory=FALSE,verbose=FALSE, ...)
{
	if(disable_cl)
	{
		nodes <- 1
	} else
	{
		if (is.null(cl)) {
			cl <- getCluster()
			on.exit( returnCluster() )
		}
	}
	
	# mask_hpc won't work if one of the two files is in memory and the other is on disk, 
	# so we will auto-write the file to disk in this case.
	
	if(xor(inMemory(x),inMemory(mask)))
	{
		filename <- tempfile()
		if(inmemory)
		{
			if(!inMemory(x))
			{
				x=readAll(x)
			} else
			{
				mask=readAll(mask)
			}
		} else
		{
			if(inMemory(x))
			{
				x=writeRaster(x,filename=filename)
			} else
			{
				mask=writeRaster(mask,filename=filename)
			}
		}
	}
	
	return(calc_hpc(stack(mask,x),
			fun=function(x)
			{ 
				nlayers_x=nlayers(x)
				mask=raster(x,layer=1)
				if(nlayers_x > 1)
				{
					pos=2:nlayers_x
					x_image=stack(mapply(function(band,inbrick)
						{raster(inbrick,layer=band)},band=pos,MoreArgs=list(inbrick=x)))
				} else
				{
					x_image=raster(x,layer=2)
				}
				x_image[mask != 1] <- NA
				return(x_image)
			},...
	,verbose=verbose))
}