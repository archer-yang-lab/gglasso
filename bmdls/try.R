

x,y,group,
				nlambda=100,lambda.min=ifelse(nobs<nvars,5e-2,1e-3),lambda, 
				standardize=TRUE,center=TRUE,eps=1e-4, 
				dfmax=as.integer(max(group))+1,pmax=min(dfmax*1.2,as.integer(max(group))),
				pf=rep(1,as.integer(max(group))),maxit=100,
				penalty.type=c("lasso","scad","mcp","sica"),
				a=ifelse(penalty.type=="scad",3.7,2.0))
{
	#################################################################################	
	#data setup
	this.call=match.call()
	if (is.null(x)) stop("Must supply x")
	if (is.null(y)) stop("Must supply y")
	np=dim(x)
  	nobs=as.integer(np[1])
  	nvars=as.integer(np[2])
  	vnames=colnames(x)
  	if(is.null(vnames))vnames=paste("V",seq(nvars),sep="")
	if(!is.numeric(y)) stop("The response must be numeric. Factors must be converted to numeric")
	if(length(y)!=nobs) stop("x and y have different number of rows")
	#################################################################################	
	#group setup
	if (!is.null(group))
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
	if(center==TRUE) 
	{
		one=rep(1, nobs)
		meanx=drop(one %*% x)/nobs
		x=scale(x, meanx, FALSE)
	}
	#################################################################################	
	#blockwise orthonormalization
	maj <- rep(0,bn)
	if(standardize==TRUE) 
	{
		maj = rep(1,bn)
		for (g in 1:bn)
		{
			ind=ix[g]:iy[g]
			decomp <- qr(x[,ind])
			if(decomp$rank < bs[g]) warning("Block belonging to columns ",  ## Warn if block has not full rank
			paste(ind, collapse = ", ")," has not full rank! \n")
			x[,ind] <- qr.Q(decomp)
		}
	}
	if(standardize==FALSE)
	{
		for(g in 1:bn) maj[g] <- max(eigen(crossprod(x[,ix[g]:iy[g]]))$values)
	}
	maj=as.double(maj)
	#################################################################################	
	#parameter setup
	if(length(pf)!=bn) stop("The size of penalty factor must be same as the number of groups")
	penalty.type=match.arg(penalty.type)
	if(penalty.type=="scad" && a<=2) stop("a must be greater than 2 for SCAD penalty")
	maxit=as.integer(maxit)
	pf=as.double(pf)
	eps=as.double(eps)
	dfmax=as.integer(dfmax)
	pmax=as.integer(pmax)
	a=as.double(a)
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
	if(penalty.type=="lasso") fit=.Fortran("blslasso",bn,bs,ix,iy,maj,
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
	if(penalty.type=="scad") fit=.Fortran("blsscad",a,bn,bs,ix,iy,maj,
								nobs,nvars,as.double(x),as.double(y),
								dfmax,pmax,nlam,flmin,ulam,eps,maxit,
								nalam=integer(1),
								b0=double(nlam),
								beta=double(nvars*nlam),
								idx=integer(pmax),
								nbeta=integer(nlam),
								alam=double(nlam),
								npass=integer(1),
								jerr=integer(1))
	if(penalty.type=="mcp") fit=.Fortran("blsmcp",a,bn,bs,ix,iy,maj,
								nobs,nvars,as.double(x),as.double(y),
								dfmax,pmax,nlam,flmin,ulam,eps,maxit,
								nalam=integer(1),
								b0=double(nlam),
								beta=double(nvars*nlam),
								idx=integer(pmax),
								nbeta=integer(nlam),
								alam=double(nlam),
								npass=integer(1),
								jerr=integer(1))
	if(penalty.type=="sica") fit=.Fortran("blssica",a,bn,bs,ix,iy,maj,
								nobs,nvars,as.double(x),as.double(y),
								dfmax,pmax,nlam,flmin,ulam,eps,maxit,
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
	print(fit)
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
		ja=fit$idx[seq(nbetamax)]#confusing but too hard to change
		oja=order(ja) 
		nzwhich=vector()
		for(j in 1:nbetamax)
		{
			g=ja[oja][j]
			nzwhich=c(nzwhich,c(ix[g]:iy[g]))
		}
		beta=matrix(fit$beta[seq(nvars*nalam)],nvars,nalam)[nzwhich,,drop=FALSE]
		df=apply(abs(beta)>0,2,sum)
		nzmax=max(df)
		ja=rep(nzwhich,nalam)
		idx=cumsum(c(1,rep(nzmax,nalam)))
		beta=new("dgCMatrix",Dim=dd,Dimnames=list(vnames,stepnames),x=as.vector(beta),p=as.integer(idx-1),i=as.integer(ja-1))
	}
	else 
	{
		beta = zeromat(nvars,nalam,vnames,stepnames)
		df=rep(0,nalam)
	}
	b0=fit$b0[seq(nalam)]
	names(b0)=stepnames
	outlist=list(b0=b0,beta=beta,df=df,
	lambda=lam,npasses=fit$npass,jerr=fit$jerr,dim=dd,call=this.call)
	class(outlist)=c("bmd.ls","bmd")
	outlist
}