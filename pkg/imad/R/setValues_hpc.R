#' @export

setValues_hpc <- function(x,values,layer=1,verbose=FALSE)
{
	require("raster")
#	nrow_x=nrow(x)
#	ncol_x=ncol(x)
#	nlayers_x=nlayers(x)
	
	if(class(values)=="numeric" || class(values)=="logical")
	{
		if(length(values)==ncell(x))
		{
		if(verbose) { print("setValues with vector...") }
		#Same as setValues
			if(missing(layer))
			{
				out=setValues(x,values)
			} else
			{
				out=setValues(x,values,layer)
			}
		} else
		{
			if(verbose) { print("setValues with longer vector...") }
			if(length(values)!=ncell(x)*nlayers(x))
			{
				stop ("values of type numeric must have a length of ncell(x)*nlayers(x)")
			} else
			{
				startidx=((c(1:nlayers(x))-1)*ncell(x))+1
				endidx=(c(2:(nlayers(x)+1))-1)*ncell(x)
				
				if(verbose) { print(startidx) }
				if(verbose) { print(endidx) }
				if(verbose) { print(length(values)) }
				
				for(i in 1:nlayers(x))
				{
					out <- setValues(x,values[startidx[i]:endidx[i]],layer=i)
				}
			}
			
		}
	}
	
	if(class(values)=="array")
	{
		if(verbose) { print("setValues with array...") }
		array_dims=dim(values)
		x_dims=c(nrow(x),ncol(x),nlayers(x))
		if(length(array_dims) !=3 || sum(array_dims==x_dims) !=3)
		{
			stop ("values of type array must have three dimensions of c(nrow(x),ncol(x),nlayers(x))")
		}
		# For-next loop not the most elegant...
		for(i in 1:nlayers(x))
		{
			out <- setValues(x,as.vector(values[,,i]),layer=i)			
		}	
		
	}
	
	if(class(values)=="matrix")
	{
		if(verbose) { print("setValues with matrix...") }
		matrix_dims=dim(values)
		x_dims=c(nrow(x),ncol(x))
		if(length(matrix_dims) !=2 || sum(matrix_dims==x_dims) !=2)
		{
			stop ("values of type matrix must have two dimensions of c(nrow(x),ncol(x))")
		}
		out <- setValues(x,as.vector(values),layer=layer)
	}
	return(out)
}