#' iMad (optimized port)
#' 
#' Perform iteratively re-weighted multivariate alteration detection.
#' 
#' @param inDataSet1 A Raster* object of the first image.
#' @param inDataSet2 A Raster* object of the second image.
#' @param pos Integer vector.  A vector of bands to use from each image.  Default is use all bands.
#' @param mask1 (Optional) A Raster* object representing a mask to be used for inDataSet1 or a numeric value to be used as the mask value.
#' @param mask2 (Optional) A Raster* object representing a mask to be used for inDataSet2 or a numeric value to be used as the mask value.
#' @param mask1_band (Optional) The band from inDataSet1 to use for masking (only if class(mask1)=="numeric").
#' @param mask2_band (Optional) The band from inDataSet2 to use for masking (only if class(mask1)=="numeric").
#' @param output_basename Character. The basename (including path) for the output files.
#' @param format Character. The output format of the rasters (see ?writeFormats).  Default is "raster".
#' @param maxiter Numeric (>= 1). The maximum number of iterations.  Default is 100.
#' @param lam Numeric. The penalization function.  CURRENTLY UNSUPPORTED.
#' @param corr_thresh Numeric. Used for situations where the canonical correlates are all nearly 1.0 (how close to 1.0 does it need to be to stop). 
#' @param delta Numeric. The smallest change in canonical correlates to end the program.
#' @param auto_extract_overlap Logical. Extract the overlap zones between the images?
#' @param reuse_existing_raster Logical. If the algorithm detects pre-create overlaps, use them?
#' @param force_extent Logical. Attempt to force the input files (and masks) to be the same extent?
#' @param enable_snow Logical. Use clusterR to (potentially) speed up calculations on a cluster? Default=FALSE.  EXPERIMENTAL.
#' @param cl Cluster.  If not assigned, the program will attempt to figure it out.  Use beginCluster() to create a cluster.
#' @param verbose Logical. Print out debugging information?
#' @param ... Passed to various raster functions (see writeRaster). Important ones include format= and overwrite=TRUE/FALSE.  datatype should be left as 'FLT4S' for proper functioning.
#' @return Returns a RasterStack object where the first layer is the chisquare image, and the subsequent layers are the iMad layers.
#' @author Mort Canty (original code) and Jonathan A. Greenberg (R port).
#' @seealso \code{\link[raster]{writeRaster}}
# @keywords 
#' @references
#' \itemize{
#' \item {Canty, M.J. and A.A. Nielsen. 2008. Automatic radiometric normalization of multitemporal satellite imagery with the iteratively re-weighted MAD transformation. Remote Sensing of Environment 112:1025-1036.}
#' \item {Nielsen, A.A. 2007. The regularized iteratively reweighted MAD method for change detection in multi- and hyperspectral data. IEEE Transactions on Image Processing 16(2):463-478.}
#' }
# @examples
# 
# \dontrun{
# } 
#' @export

# TODO: parallelize if/where possible
# TODO: integrate penalization function for hyperspectral imagery.
# TODO: write out full mask and check for existing masks.
# TODO: give the output stack names
# TODO: use layerStats instead of cov.wt.raster
# TODO: allow mask to be a single value
# TODO: subset bands
# TODO: save the mask step.

iMad <- function(inDataSet1,inDataSet2,pos,
		mask1,mask2,mask1_band=1,mask2_band=1,
		output_basename, format="raster",
		maxiter=100,lam=0,delta=0.001,corr_thresh=0.001,
		auto_extract_overlap=TRUE,reuse_existing_raster=TRUE,force_extent=TRUE,
		enable_snow=FALSE,cl=NULL,
		verbose=FALSE,
		...)
{
	require("raster")		
	if(verbose) { print("Verbose mode enabled...")}
	
	# Do some pre-checks up here.
	if(verbose) { print("Initial checks...")}
	if(!inherits(inDataSet1,'Raster')) stop ("inDataSet1 must be a Raster* object...")
	if(!inherits(inDataSet2,'Raster')) stop ("inDataSet2 must be a Raster* object...")
	if(maxiter < 1) stop ("maxiter must be greater than 0...")
	if(missing(output_basename)) stop ("output_basename must be set...")
	
	# End pre-checks.

	# Setup output filenames.
	output_inDataSet1_subset_filename=paste(output_basename,"_inDataSet1_overlap",sep="")
	output_inDataSet2_subset_filename=paste(output_basename,"_inDataSet2_overlap",sep="")
	output_inDataSet1_subset_mask_filename=paste(output_basename,"_inDataSet1_mask_overlap",sep="")
	output_inDataSet2_subset_mask_filename=paste(output_basename,"_inDataSet2_mask_overlap",sep="")
	output_inDataSet1_masked=paste(output_basename,"_inDataSet1_masked",sep="")
	output_inDataSet2_masked=paste(output_basename,"_inDataSet2_masked",sep="")
	
	
	output_MAD_filename=paste(output_basename,"_iMAD",sep="")
	output_chisqr_filename=paste(output_basename,"_iMAD_chisqr",sep="")
	
##### Subset out bands if needed.
	if(!missing(pos))
	{
		if(verbose) { print("Subsetting bands...")}
#		inDataSet1=stack(inDataSet1,bands=pos)
#		inDataSet2=stack(inDataSet2,bands=pos)
		inDataSet1=stack(mapply(function(band,inbrick){raster(inbrick,layer=band)},band=pos,MoreArgs=list(inbrick=inDataSet1)))
		inDataSet2=stack(mapply(function(band,inbrick){raster(inbrick,layer=band)},band=pos,MoreArgs=list(inbrick=inDataSet2)))
	}
	
	bands=nlayers(inDataSet1)
	pos=0:(bands-1)
	
##### Extract overlap regions
	if(auto_extract_overlap)
	{
		if(verbose) { print("Extracting the overlap region...") }

		# Not working yet, needs to check for the the filename + extension.
		
		if(missing(format))
		{
			format="raster"
		}
		
		output_inDataSet1_subset_full_filetype <- raster:::.filetype(format=format, filename=output_inDataSet1_subset_filename)
		output_inDataSet2_subset_full_filetype <- raster:::.filetype(format=format, filename=output_inDataSet2_subset_filename)
		
		output_inDataSet1_subset_full_filename <- raster:::.getExtension(output_inDataSet1_subset_filename, output_inDataSet1_subset_full_filetype)
		output_inDataSet2_subset_full_filename <- raster:::.getExtension(output_inDataSet2_subset_filename, output_inDataSet2_subset_full_filetype)
		
		if(!reuse_existing_raster || 
				!file.exists(output_inDataSet1_subset_full_filename) || !file.exists(output_inDataSet2_subset_full_filename))
		{
		# If the user does not want to reuse existing cropped datasets or if they don't exist...
			# First check to see if they line up...
		
			# If not, run the extraction...
#			overlaps=extract_overlap_rasters(inDataSet1,inDataSet2,
#					filename1=output_inDataSet1_subset_filename,filename2=output_inDataSet2_subset_filename,
#					verbose=verbose,format=format,datatype='FLT4S',...)
#			inDataSet1=overlaps[[1]]
#			inDataSet2=overlaps[[2]]
			inDataSet1=crop(inDataSet1,inDataSet2,filename=output_inDataSet1_subset_filename,datatype='FLT4S',verbose=verbose,...)
			inDataSet2=crop(inDataSet2,inDataSet1,filename=output_inDataSet2_subset_filename,datatype='FLT4S',verbose=verbose,...)
		} else
		{
		# Otherwise use the existing cropped datasets.
			print("Reusing the existing regions...")
			inDataSet1=brick(output_inDataSet1_subset_full_filename)
			inDataSet2=brick(output_inDataSet2_subset_full_filename)
		}		
		
		# Now crop the masks and add them together.
		if((!missing(mask1) || !missing(mask2)) && auto_extract_overlap)
		{
#			if(verbose) { print(class(mask1)) }
			if(!missing(mask1) && !(class(mask1)=="numeric"))
			{
				if(verbose) { print("Cropping mask1...") }
				mask1=crop(mask1,inDataSet1)
#				mask1_overlap=extract_overlap_rasters(mask1,inDataSet2,
#						filename1=output_inDataSet1_subset_mask_filename,
#						raster2_crop=FALSE,
#						verbose=verbose,format=format,...)
#				mask1=mask1_overlap
			}
			
			if(!missing(mask2) && !(class(mask2)=="numeric"))
			{
				if(verbose) { print("Cropping mask2...") }
				mask2=crop(mask2,inDataSet2)
#				mask2_overlap=extract_overlap_rasters(inDataSet1,mask2,
#						filename2=output_inDataSet2_subset_mask_filename,
#						raster1_crop=FALSE,
#						verbose=verbose,format=format,...)
#				mask2=mask2_overlap
			}
		}
		
		
	} else
	{
		if(verbose) { print("Not extracting the overlap region...") }
		# Should check for overlap here...
	}

	
##### Fix extents if need be...

	# Fix for extents that differ.  We should combine this with a col/row check.
	if(force_extent)
	{
		if(verbose) { print("Forcing extents...") }
		extent(inDataSet2)=extent(inDataSet1)
	}
	
##### Set up the masks
	mask=NA
	# Check for numeric masks
	if(!missing(mask1))
	{
		if(class(mask1)=="numeric")
		{
			if(verbose) { print("Creating mask1...") }
			### TODO: HPC
			if(enable_snow)
			{
				mask1=calc_hpc(x=raster(inDataSet1,layer=mask1_band),function(x,mask1) { x != mask1 },args=list(mask1=mask1))			
			} else
			{
				mask1=raster(inDataSet1,layer=mask1_band) != mask1
			}
		}
	}
	
	if(!missing(mask2))
	{
		if(class(mask2)=="numeric")
		{
			if(verbose) { print("Creating mask2...") }
			if(enable_snow)
			{
				mask2=calc_hpc(x=raster(inDataSet2,layer=mask2_band),function(x,mask2) { x != mask2 },args=list(mask2=mask2))
			} else
			{
				mask2=raster(inDataSet2,layer=mask2_band) != mask2
			}
		}
	}
	
	if(!missing(mask1) && !missing(mask2))
	{
		if(verbose) { print("Both masks present...") }
		if(force_extent)
		{
			extent(mask1)=extent(inDataSet1)
			extent(mask2)=extent(inDataSet1)
		}
		# This can be HPCed
		mask=mask1*mask2
	} else
	{
		if(!missing(mask1))
		{
			if(verbose) { print("Mask #1 present...") }
			if(force_extent)
			{
				extent(mask1)=extent(inDataSet1)
			}
			mask=mask1
		} else
		{
			if(verbose) { print("Mask #2 present...") }
			if(force_extent)
			{
				extent(mask2)=extent(inDataSet1)
			}
			mask=mask2
		}
		# This can be HPCed
		mask[mask==0] <- NA
	}

	cols=ncol(inDataSet1)
	rows=nrow(inDataSet1)
	
	cols2=ncol(inDataSet2)
	rows2=nrow(inDataSet2)
	
	if(cols != cols2 || rows != rows2) stop("Input rows and columns must be the same, try using auto_extract_overlap=TRUE...")
	
	if(class(mask)!="logical")
	{
		if(verbose) { print("Masking...") }
		if(file.exists(build_raster_filename(output_inDataSet1_masked,format=format)) && reuse_existing_raster)
		{
			inDataSet1=brick(build_raster_filename(output_inDataSet1_masked,format=format))
		} else
		{
			if(enable_snow)
			{
				if(verbose) { print("HPC masking inDataSet1...") }
				inDataSet1=mask_hpc(inDataSet1,mask)
			} else
			{
				inDataSet1=mask(x=inDataSet1,mask=mask,filename=output_inDataSet1_masked,...)
			}
		}

		if(file.exists(build_raster_filename(output_inDataSet2_masked,format=format)) && reuse_existing_raster)
		{
			inDataSet2=brick(build_raster_filename(output_inDataSet2_masked,format=format))
		} else
		{
			if(enable_snow)
			{
				if(verbose) { print("HPC masking inDataSet2...") }
				inDataSet2=mask_hpc(inDataSet2,mask)
			} else
			{
				inDataSet2=mask(x=inDataSet2,mask=mask,filename=output_inDataSet2_masked,...)
			}
		}

#		inDataSet1 <- inDataSet1*mask
#		inDataSet2 <- inDataSet2*mask
	}
	
	if(verbose) { print("Creating initial weighting raster and stacking inDataSets...")}
	if(enable_snow)
	{
		wt = calc_hpc(raster(inDataSet1,layer=1),fun=function(x) { x*0+1 })	
		
		wt = mask_hpc(wt,mask)
	} else
	{
		wt = raster(inDataSet1,layer=1)*0+1
	}
	
	dm = stack(inDataSet1,inDataSet2)
	
	delta = 1.0
	oldrho = array(data=0,dim=bands)
	iter = 0
	ab_nan=FALSE

# Mods to include the penalization function.  Comment this out if this chokes.
#	if(lam>0) { Omega_L = diag(bands) }

	
### MAIN LOOP
	while(delta > 0.001 && iter < maxiter && !(ab_nan))
	{
		if(verbose)
		{
			print(paste("Iteration:",iter))
		}
		
		# This needs to be swapped with "layerStats" in raster
		
#		sigma_means=cov.wt.raster(dm,wt)

		if(enable_snow)
		{
			if(verbose) { print("HPC calculating weighted covariance and means...")}
			sigma_means=layerStats_hpc(x=dm,w=wt,stat='weighted.cov',na.rm=TRUE,enable_snow=TRUE,verbose=verbose)
		} else
		{
			if(verbose) { print("Calculating weighted covariance and means...")}
			sigma_means=layerStats(dm,'weighted.cov',wt,na.rm=TRUE)
		}	
		
		if(verbose) { print(sigma_means) }
		
		sigma=sigma_means[[1]]
		means=sigma_means[[2]]
		
		if(is.nan(mean(sigma)) || is.nan(mean(means)))
		{
			if(verbose)
			{
				print("Possible numerical precision problem with sigma or means, exiting loop and using the current results...")
				ab_nan=TRUE
			}
		} else
		{	
			s11 = sigma[1:(bands),1:(bands)]
			s22 = sigma[(bands+1):(2*bands),(bands+1):(2*bands)]
			s12 = sigma[1:(bands),(bands+1):(2*bands)]
			s21 = sigma[(bands+1):(2*bands),1:(bands)]
			
			# Mods to include the penalization function.  Comment this out if this chokes.
		#		s11 = (1-lam)*S11 + lam*Omega_L
		#		s22 = (1-lam)*S22 + lam*Omega_L
#			if(lam>0) {
#				s11 = (1-lam)*S11 + lam*Omega_L
#				s22 = (1-lam)*S22 + lam*Omega_L
#			}

			if(verbose) { print("Calculating generalized eigenvalues and eigenvectors...")}		
			lama_a=Rdggev(JOBVL=F,JOBVR=T,A=s12%*%solve(s22)%*%s21,B=s11)
			a=lama_a$VR
				lama=lama_a$GENEIGENVALUES
				
			lamb_b=Rdggev(JOBVL=F,JOBVR=T,A=s21%*%solve(s11)%*%s12,B=s22)
			b=lamb_b$VR
			lamb=lamb_b$GENEIGENVALUES
	
			# Eigenvalues, ranked largest to smallest
			idx=rank(lama)
			a=a[,idx]
				
			idx=rank(lamb)
			b=b[,idx]
	
	# Penalization stuff needs to go somewhere around here
	# IDL Code:
	
	
	#
	#; stopping criterion
	#delta = max(abs(rho-old_rho))
	#old_rho = rho
	#
	#; ensure sum of positive correlations between X and U is positive
	#; their covariance matrix is S11##A
	#invsig_x = diag_matrix(1/sqrt(diag_matrix(S11)))
	#sum = total(invsig_x##S11##A,2)
	#				A = A##diag_matrix(sum/abs(sum))   
	#
	#; ensure positive correlation between each pair of canonical variates
	#cov = diag_matrix(transpose(A)##S12##B)
	#B = B##diag_matrix(cov/abs(cov))
	#
	#if iter gt 0 and iter eq niter then goto, done 
			
			#mu = sqrt(mu2)  
			rho=sqrt(lamb[idx])	
#			mu=rho
			
			# Fix for near-perfect correlations
			if(abs(sum(rho)-bands) < corr_thresh*bands)
			{
				print("Perfect correlation, exiting...")
				ab_nan=TRUE
			} else
			{
				# normalize dispersions   
	
	#a2=diag_matrix(transpose(A)##A)
#				a2=diag(t(a)%*%a)
#				b2=diag(t(b)%*%b)
#				sigma2=sqrt((2-lam*(a2*b2))/(1-lam-2*mu))
#				rho=mu*(1-lam)/sqrt((1-lam*a2)*(1-lam*b2))
#				sigMads
	#b2=diag_matrix(transpose(B)##B) 
	#sigma = sqrt( (2-lam*(a2+b2))/(1-lam)-2*mu )
	#rho=mu*(1-lam)/sqrt( (1-lam*a2)*(1-lam*b2) )     
	#
	#sigMads = ones##sigma 
	#means1  = ones##means[0:num_bands-1]
	#means2  = ones##means[num_bands:*]
	
				tmp1=t(a)%*%s11%*%a
				tmp2=1/(sqrt(diag(tmp1)))
				tmp3=t(array(tmp2,c(bands,length(tmp2))))
				a=a*tmp3
					
				tmp1=t(b)%*%s22%*%b
				tmp2=1/(sqrt(diag(tmp1)))
				tmp3=t(array(tmp2,c(bands,length(tmp2))))
				b=b*tmp3
		
				# assure positive correlation
				tmp=diag(t(a)%*%s12%*%b)
				b=b%*%diag(tmp/abs(tmp))
							
				if(is.nan(mean(a)) || is.nan(mean(b)))
				{
					if(verbose)
					{
						print("Possible numerical precision problem with a or b, exiting loop and using the current results...")
						ab_nan=TRUE
					}
				} else
				{
					#     canonical and MAD variates
					if(verbose) { print("Calculating canonical and MAD variates...")}
					means_a=means[1:bands]
					# Experimental clusterR
	
					if(verbose) { print(a) }
					if(verbose) { print(means_a) }
#					if(enable_snow)
#					{
##						U=calc_hpc(x=inDataSet1,fun=function(x,a,means_a) { as.vector(t(a)%*%(x-means_a)) }, 
##								args=list(a=a,means_a=means_a))
#						
#						U=calc_hpc(x=inDataSet1,fun=function(x,a,means_a) { as.vector(a*(x-means_a)) }, 
#								args=list(a=a,means_a=means_a))
#						
##						U=clusterR(x=inDataSet1,fun=function(x,a,means_a) { as.vector(t(a)%*%(x-means_a)) }, 
##								args=list(a=a,means_a=means_a) )
#					} else
#					{
						U=calc(inDataSet1,fun=function(x) { as.vector(t(a)*(x-means_a)) } )
#					}
					
					means_b=means[(bands+1):(bands*2)]
#					if(enable_snow)
#					{
##						V=calc_hpc(x=inDataSet2,fun=function(x,b,means_b) { as.vector(t(b)%*%(x-means_b)) }, 
##								args=list(b=b,means_b=means_b) )
#						
#						V=calc_hpc(x=inDataSet2,fun=function(x,b,means_b) { as.vector(b*(x-means_b)) }, 
#								args=list(b=b,means_b=means_b) )
#					} else
#					{
						V=calc(inDataSet2,fun=function(x) { as.vector(t(b)%*%(x-means_b)) } )
#					}
					
#					if(enable_snow)
#					{
#						MAD=calc_hpc(x=stack(U,V),
#							fun=function(x,b1,b2)
#							{
#								(raster(x,layer=b1)-raster(x,layer=b2))
#							},args=list(b1=1,b2=2))
#					} else
#					{
						MAD = U-V
#					}	
					#     new weights	
					if(verbose) { print("Generating new weights...")}
					var_mad=t(2*(1-rho))
					chisqr=calc(MAD,fun=function(x) { sum(x^2/var_mad) })
					wt=1-calc(chisqr,fun=function(x) { pchisq(x,bands) })
						
					delta = sum(abs(rho-oldrho))
					oldrho = rho
					if(verbose)
					{
						print(paste("Delta:",delta)) 
						print("rho:")
						print(rho)
						print("****************")
					}
					iter = iter+1
				}
			}
		}
	}
	
	# Output results
	MAD_brick=writeRaster(MAD,filename=output_MAD_filename,format=format,...)
	chisqr_raster=writeRaster(chisqr,filename=output_chisqr_filename,format=format,...)
	
	return(stack(chisqr_raster,MAD_brick))
	
}