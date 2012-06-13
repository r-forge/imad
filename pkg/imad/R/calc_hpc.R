# clusterR using mmap
# Original code by Robert Hijimans, mmap integration by Jonathan Greenberg
#' @export

calc_hpc <- function(x, fun, args=NULL, filename='', cl=NULL, m=2, disable_cl=FALSE,
		todisk=FALSE,verbose=FALSE,...) 
{
	require("raster")
	require("snowfall")
		
	# Do some file checks up here.
	
	if(verbose) { print("Setting up cluster...")}
	
	if(disable_cl)
	{
		nodes <- 1
	} else
	{
		if (is.null(cl)) {
			# Check to see if a cluster is running
#			if(sfIsRunning())
#			{
#				cl <- getCluster()
#				cluster_shutdown=FALSE
#			} else
#			{
#				cl <- beginCluster()
				cl <- getCluster()
#				cluster_shutdown=TRUE
#			}
		}
		nodes <- length(cl)
	}

	# We should test a single pixel here to see the size of the output...
	
	# We are going to pull out the first row and first two pixels to check the function...
	if(verbose) { print("Performing a pre-check...")}
	r_check <- crop(x, extent(x, r1=1, r2=1, c1=1,c2=2))
	
	if(!is.null(args)) {
		r_check_function <- do.call(fun, c(r_check, args))
		if(inherits(r_check_function,"Raster")) { r_check_function = getValues(r_check_function) }
	} else
	{
		r_check_function <- fun(r_check)
		if(inherits(r_check_function,"Raster")) { r_check_function = getValues(r_check_function) }
	}
	
	# Next we do a check for the number of output bands.  A matrix output indicates a multi-band output.
	if(verbose) { print("Determining number of output bands...")}
	if(class(r_check_function)=="numeric" || class(r_check_function)=="logical")
	{
		outbands=1
	} else
	{
		outbands=dim(r_check_function)[2]
	}
	if(verbose) { print(outbands) }
	
	# The algorithm works differently if it processes in memory, so we setup the output here.	
	if(canProcessInMemory(raster(x),n=outbands) && !todisk)
	{
		if(verbose) { print("Processing in memory...") }
		inmemory=TRUE
		if(outbands > 1)
		{
			outraster <- brick(x,nl=outbands)
#			outraster <- readAll(outraster)
		} else
		{
			outraster <- raster(x)
#			outraster <- readAll(outraster)
		}
#		out<-array(dim=c(nrow(x),ncol(x),outbands))
	} else
	{
		if(verbose) { print("Not processing in memory...") }
		inmemory=FALSE
		if(verbose) { print("Creating output file with ff...")}
		require("ff")
		require("mmap")
		
		outdata_ncells=nrow(x)*ncol(x)*outbands
		if(verbose) { print(outdata_ncells) }
		
		if(filename=="")
		{	
			filename <- tempfile()
		} 
		if(verbose) { print(filename) }
		# How about using ff?
		out<-ff(vmode="double",length=outdata_ncells,filename=filename)
		finalizer(out) <- close
		close(out)	
	}
	
	if(verbose) { print("Determining optimal block size...")}
	m <- max(1, round(m))
	tr <- blockSize(x, minblocks=nodes*m )
	if (tr$n < nodes) 
	{
		nodes <- tr$n
	}
	
	tr$row2 <- tr$row + tr$nrows - 1

	tr_out=list(row=((tr$row-1)*outbands+1))
	tr_out$row2=((tr$row2)*outbands)
	
	i=1:tr$n
	
#	if(disable_cl)
#	# Use only for debugging.
#	{
#		mapply(function(i,fun,args,x,tr,filename,outbands) 
#			{
#				r <- crop(x, extent(x, r1=tr$row[i], r2=tr$row2[i], c1=1, c2=ncol(x)))
#				if(!is.null(args)) {
#					r <- fun(r) 
#				} else
#				{
#					r <- do.call(fun, c(r, args))
#				}
#				out <- mmap(filename,mode=real64())
#				cellStart=((cellFromRowCol(x,row=tr$row[i],col=1))-1)*outbands+1
#				cellEnd=((cellFromRowCol(x,row=tr$row2[i],col=ncol(x))))*outbands
#				out[cellStart:cellEnd] <- as.vector(t(getValues(r)))
##				out[cellFromRow(x,tr$row[i]:tr$row2[i])] <- as.vector(getValues(r))
#				munmap(out)
#				return(NULL)
#			},
#			i,MoreArgs=list(fun=fun,x=x,tr=tr,args=args,filename=filename,outbands=outbands))
#
#	} else
#	{
		if(verbose) { print("Starting the cluster function...")}
		out <- clusterMap(cl,function(fun,i,args,x,tr,filename,outbands,inmemory,verbose) 
			{
				r <- crop(x, extent(x, r1=tr$row[i], r2=tr$row2[i], c1=1, c2=ncol(x)))
				if(is.null(args)) {
					r <- fun(r) 
				} else
				{
					r <- do.call(fun, c(r, args))
				}
#				if(verbose) { print(class(r)) }
#				if(verbose) { print(dim(r)) }
				
				if(inmemory)
				{
					if(inherits(r, 'Raster'))
					{
						r=getValues(r)
					}
					return(r)
				} else
				{
					# This performs parallel writes
					cellStart=((cellFromRowCol(x,row=tr$row[i],col=1))-1)*outbands+1
					cellEnd=((cellFromRowCol(x,row=tr$row2[i],col=ncol(x))))*outbands
					# Disable transpose for BIL?
					out <- mmap(filename,mode=real64())
					if(inherits(r, 'Raster')) 
						{ out[cellStart:cellEnd] <- as.vector(t(getValues(r))) }
					else
						{ out[cellStart:cellEnd] <- as.vector(t(r)) }
					munmap(out)
					return(NULL)
				}
			},
			i,MoreArgs=list(fun=fun,x=x,tr=tr,args=args,filename=filename,outbands=outbands,inmemory=inmemory,
					verbose=verbose))
#	}
		
	if(inmemory)
	{
#		if(verbose) { print(length(out[[1]])) }
#		print(class(unlist(out)))
#		print(dims(unlist(out)))
		outraster=setValues_hpc(outraster,unlist(out),verbose=verbose)
#		print(outraster_test)
	} else
	{
		# Let's see if we can trick raster into making us a proper header...
		if(outbands > 1) 
		{ 
			reference_raster=brick(raster(x,layer=1),nl=outbands) 
		} else
		{
			if(nlayers(x) > 1) { reference_raster=raster(x,layer=1) } else
			{ reference_raster=x }	
		}
		outraster_base <- writeStart(reference_raster,filename=paste(filename,".grd",sep=""),datatype="FLT8S",bandorder="BIP",...)
		suppressWarnings(outraster_base <- writeStop(outraster_base))
		file.remove(paste(filename,".gri",sep=""))
		file.rename(filename,paste(filename,".gri",sep=""))
		outraster=brick(paste(filename,".grd",sep=""))
	}
	
#	if(cluster_shutdown) { endCluster() }
	# Some housekeeping
	out <- clusterMap(cl,function(x) { x },x=0)
	
	return(outraster)
}

