#' @export

brick_to_stack <- function(inbrick)
{
	brick_nlayers=nlayers(inbrick)
	pos=1:brick_nlayers
	return(
		stack(
			mapply(function(band,inbrick){ raster(inbrick,layer=band) },band=pos,MoreArgs=list(inbrick=inbrick),SIMPLIFY=FALSE))
	)	
}