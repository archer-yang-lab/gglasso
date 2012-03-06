bmd.logit<-function(x,y,group,
				nlambda=100,lambda.min=ifelse(nobs<nvars,5e-2,1e-3),lambda, 
				standardize=TRUE,eps=1e-4, 
				dfmax=as.integer(max(group))+1,pmax=min(dfmax*1.2,as.integer(max(group))),
				pf=rep(1,as.integer(max(group))),maxit=100)
{
	#################################################################################	
	#data setup
	this.call=match.call()
	np=dim(x)
  	nobs=as.integer(np[1])
  	nvars=as.integer(np[2])
  	vnames=colnames(x)
  	if(is.null(vnames))vnames=paste("V",seq(nvars),sep="")
	if(!is.numeric(y)) stop("The response must be numeric. Factors must be converted to numeric")
	if(length(y)!=nobs) stop("x and y have different number of rows")
    if(!all(is.element(y,c(-1,1)))) stop("classification requires the response to be in {-1,1}")
	#################################################################################	
	#group setup
	if (!missing(group))
    {
        if (length(group)!=nvars) stop("group does not match x")
    }else group=1:nvars
	bn=as.integer(max(group))
    bs=as.integer(as.numeric(table(group)))
    if(!identical(as.integer(sort(unique(group))),as.integer(1:bn))) stop("Groups must be consecutively numbered 1,2,3,...")
	ix=rep(0,bn)
	iy=rep(0,bn)
	j=1
	for(g in 1:bn)
	{
		ix[g]=j
		iy[g]=j+bs[g]-1
		j=j+bs[g]
	}
	ix=as.integer(ix)
	iy=as.integer(iy)
	#################################################################################	
	#centering input variable
	one=rep(1, nobs)
	meanx=drop(one %*% x)/nobs
	x=scale(x, meanx, FALSE)
	maj <- rep(0,bn)
	if(standardize==TRUE) 
	{
		maj = rep(1,bn)
		for (g in 1:bn)
		{
			ind=ix[g]:iy[g]
			xind=as.matrix(x[,ind])
			px=ncol(xind)
			if(px>nobs) stop("Too much variables to orthogonalize for each group")
			decomp <- qr(xind)
			if(decomp$rank < bs[g]) stop("Block belonging to columns ",  ## Warn if block has not full rank
			paste(ind, collapse = ", ")," has not full rank! \n")
			x[,ind] <- qr.Q(decomp)
		}
	}
	else{
		for(g in 1:bn) maj[g] <- max(eigen(crossprod(x[,ix[g]:iy[g]]))$values)
	}
	maj=as.double(maj)
	#################################################################################	
	#parameter setup
	if(length(pf)!=bn) stop("The size of penalty factor must be same as the number of groups")
	maxit=as.integer(maxit)
	pf=as.double(pf)
	eps=as.double(eps)
	dfmax=as.integer(dfmax)
	pmax=as.integer(pmax)
	#################################################################################	
	#lambda setup
	nlam=as.integer(nlambda)
	if(missing(lambda))
	{
		if(lambda.min>=1) stop("lambda.min should be less than 1")
		if(lambda.min<1.0E-6) stop("lambda.min is too small")
		flmin=as.double(lambda.min)
		ulam=double(1) #ulam=0 if lambda is missing
	}
	else
	{
		#flmin=1 if user define lambda
		flmin=as.double(1)    
		if(any(lambda<0)) stop("lambdas should be non-negative")
		ulam=as.double(rev(sort(lambda))) #lambda is declining
		nlam=as.integer(length(lambda))
	}
	#################################################################################	
	# call Fortran core
	fit=.Fortran("log_f",bn,bs,ix,iy,maj,
								nobs,nvars,as.double(x),as.double(y),
								pf,dfmax,pmax,nlam,flmin,ulam,eps,maxit,
								nalam=integer(1),
								b0=double(nlam),
								beta=double(nvars*nlam),
								idx=integer(pmax),
								nbeta=integer(nlam),
								alam=double(nlam),
								npass=integer(1),
								jerr=integer(1))
	#################################################################################	
	# output		
	nalam=fit$nalam
	nbeta=fit$nbeta[seq(nalam)]
	nbetamax=max(nbeta)
	lam=fit$alam[seq(nalam)]
	if(missing(lambda)) lam=lamfix(lam)##first lambda is infinity; changed to entry point
	stepnames=paste("s",seq(nalam)-1,sep="")
	
	errmsg=err(fit$jerr,maxit,pmax)### error messages from fortran
	switch(paste(errmsg$n),
	"1"=stop(errmsg$msg,call.=FALSE),
	"-1"=warning(errmsg$msg,call.=FALSE)
	)
	
	dd=c(nvars,nalam)
	if(nbetamax>0)
	{
		beta=matrix(fit$beta[seq(nvars*nalam)],nvars,nalam,dimnames=list(vnames,stepnames))
		df=apply(abs(beta)>0,2,sum)
	}
	else 
	{
		beta = matrix(0,nvars,nalam,dimnames=list(vnames,stepnames))
		df=rep(0,nalam)
	}
	b0=fit$b0[seq(nalam)]
	names(b0)=stepnames
	outlist=list(b0=b0,beta=beta,df=df,
	lambda=lam,npasses=fit$npass,jerr=fit$jerr,dim=dd,call=this.call)
	class(outlist)=c("bmd.logit","bmd")
	outlist
}