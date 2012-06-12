#' iMad Original (direct python port)
#' 
#' Perform iteratively re-weighted multivariate alteration detection.
#' 
#' @param inDataSet1 A Raster* object of the first image.
#' @param inDataSet2 A Raster* object of the second image.
#' @param maxiter Numeric (>= 1). The maximum number of iterations.
#' @param lam Numeric. The penalization function.
#' @param output_basename Character. The basename (including path) for the output files.
#' @return Returns a RasterStack object where the first layer is the chisquare image, and the subsequent layers are the iMad layers.
#' @author Mort Canty (original code) and Jonathan A. Greenberg (R port).
# @seealso 
# @keywords 
#' @examples
#' 
#' \dontrun{
#' } 

iMad_original <- function(inDataSet1,inDataSet2,maxiter=100,lam=0,output_basename,verbose=FALSE,...)
{
	require("raster")	
	
	covw <- function(dm,w)
	{
		N = dim(dm)[1]
		n = dim(w)
		sumw = sum(w)
		ws = array(data=w,dim=c(N,n))
		means = array(data=rowSums(ws*dm)/sumw,c(N,1))
		means = array(data=means,c(N,n))
		dmc = dm - means
		dmc = dmc*sqrt(ws)
		covmat = dmc%*%t(dmc)/sumw
		return(list(covmat,means))
	}
	
	cols=ncol(inDataSet1)
	rows=nrow(inDataSet1)
	bands=nlayers(inDataSet1)
	pos=0:(bands-1)
	
	wt = array(data=1,dim=cols*rows)
	dm = array(data=0,dim=c(2*bands,cols*rows))
	
	k = 1
	for(b in pos)
	{
		dm[k,]=as.vector(t(as.matrix(raster(inDataSet1,layer=(b+1)))))
		dm[(bands+k),]=as.vector(t(as.matrix(raster(inDataSet2,layer=(b+1)))))
		k=k+1
	}
	
	delta = 1.0
	oldrho = array(data=0,dim=bands)
	iter = 0
#	ab_nan=FALSE
	while(delta > 0.001 && iter < maxiter)
	{
		if(verbose)
		{
			print(paste("Iteration:",iter))
		}
		
		
		#	sigma,means = covw(dm,wt)
		sigma_means=covw(dm,wt)
		
		sigma=sigma_means[[1]]
		means=sigma_means[[2]]
		
		s11 = sigma[1:(bands),1:(bands)]
		s22 = sigma[(bands+1):(2*bands),(bands+1):(2*bands)]
		s12 = sigma[1:(bands),(bands+1):(2*bands)]
		s21 = sigma[(bands+1):(2*bands),1:(bands)]
		
		lama_a=Rdggev(JOBVL=F,JOBVR=T,A=s12%*%solve(s22)%*%s21,B=s11)
		a=lama_a$VR
		lama=lama_a$GENEIGENVALUES
		
		lamb_b=Rdggev(JOBVL=F,JOBVR=T,A=s21%*%solve(s11)%*%s12,B=s22)
		b=lamb_b$VR
		lamb=lamb_b$GENEIGENVALUES
		
		idx=rank(lama)
		a=a[,idx]
		
		idx=rank(lamb)
		b=b[,idx]
		
		rho=sqrt(lamb[idx])	
		
		# normalize dispersions   
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
		
		U=t(a)%*%(dm[1:bands,]-means[1:bands,])
		V=t(b)%*%(dm[(bands+1):(bands*2),]-means[(bands+1):(bands*2),])
		MAD = U-V
			
			#     new weights
		var_mad=array(t(2*(1-rho)),c(length(rho),rows*cols))
		chisqr = colSums((MAD*MAD)/var_mad)
		wt=array(data=(1-pchisq(chisqr,bands)),dim=cols*rows)
			
		delta = sum(abs(rho-oldrho))
		oldrho = rho
		if(verbose)
		{
			print(paste("Delta:",delta)) 
			print(rho)
			print("****************")
		}
			iter = iter+1
	}
	
	# Output results
	
	output_MAD_filename=paste(output_basename,"_iMAD")
	output_chisqr_filename=paste(output_basename,"_iMAD_chisqr")
	
	MAD_image=array(MAD,c(rows,cols,bands))
	MAD_brick=brick(MAD_image)
	MAD_brick=writeRaster(MAD_brick,filename=output_MAD_filename,...)
	
	chisqr_image=matrix(data=chisqr,nrow=rows,ncol=cols)
	chisqr_raster=raster(chisqr_image)
	chisqr_raster=writeRaster(chisqr_raster,filename=output_chisqr_filename,...)
	
	return(stack(chisqr_raster,MAD_brick))
	
}