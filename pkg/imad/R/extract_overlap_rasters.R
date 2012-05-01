#' Coerces two rasters to have a shared subset of their common overlap.
#' @title extract_overlap_rasters
#' @param raster1 The first raster to extract the overlap zone from.
#' @param raster2 The second raster to extract the overlap zone from.
#' @param filename1 The output filename of raster1 (optional).
#' @param filename2 The output filename of raster2 (optional).
#' @param raster1_crop Logical. Proceed with cropping and returning raster1?
#' @param raster2_crop Logical. Proceed with cropping and returning raster2?
#' @param verbose Logical. Provide debugging feedback?
#' @param ... Other parameters to pass to crop.
#' @name extract_overlap_rasters
#' @author Jonathan A. Greenberg \email{STARStools@@estarcion.net}
#' @seealso \code{\link[raster]{crop}}
#' @export


extract_overlap_rasters <- function(raster1,raster2,filename1,filename2,raster1_crop=TRUE,raster2_crop=TRUE,verbose=TRUE,...)
{
	if(!raster1_crop && !raster2_crop) stop("raster1_crop and/or raster2_crop must be set to true")
	
	extent1=bbox(raster1)
	extent2=bbox(raster2)
	
	overlap_extent=extent1
	for(i in 1:2)
	{
	overlap_extent[i]=max(c(extent1[i],extent2[i]))
	}
	for(i in 3:4)
	{
		overlap_extent[i]=min(c(extent1[i],extent2[i]))
	}
	
	
	if(raster1_crop)
	{
		if(verbose) { print("Cropping raster1...") }
		raster1_overlap=crop(raster1,overlap_extent,filename=filename1,...)
		extent(raster1_overlap)=overlap_extent
	}
	if(raster2_crop)
	{
		if(verbose) { print("Cropping raster2...") }
		raster2_overlap=crop(raster2,overlap_extent,filename=filename2,...)
		extent(raster2_overlap)=overlap_extent
	}
	
	if(raster1_crop && raster2_crop)
	{
		return(list(raster1_overlap,raster2_overlap))
	} else
	{
		if(raster1_crop) { return(raster1_overlap) } else { return(raster2_overlap) }	
	}
}

