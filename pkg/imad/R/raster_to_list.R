#' @export

raster_to_list=function(inraster)
{
	inraster_nlayers=nlayers(inraster)
	pos=1:inraster_nlayers
	return(
		mapply(function(band,inraster){ raster(inraster,layer=band) },band=pos,MoreArgs=list(inraster=inraster),SIMPLIFY=FALSE)
	)
	
}