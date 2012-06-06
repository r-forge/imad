#' @export

spectral_subset <- function(x,pos)
{
	if(!inherits(x, 'Raster'))
	{
		stop ("x must be a Raster* object...")
	}
	
	if(class(pos)!="integer")
	{
		stop ("pos must be an integer vector...")
	}
	
	nlayers_x=nlayers(x)
	posmin=min(pos)
	posmax=max(pos)
	
	if(posmin < 1 || posmax > nlayers_x)
	{
		stop ("pos must range between 1 and nlayers(x)...")
	}
	
	x_sub=
		stack(
		mapply(
			band=pos,MoreArgs=list(inbrick=x),
			function(band,inbrick)
				{raster(inbrick,layer=band)}
		)
	)

	return(x_sub)
}