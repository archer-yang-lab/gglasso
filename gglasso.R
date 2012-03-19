gglasso <-function(x,y,group=NULL,
	           loss = c("ls","hreg","logit","sqsvm","hsvm"),
				nlambda=100,lambda.factor=ifelse(nobs<nvars,5e-2,1e-3),lambda=NULL, 
				weights,
				pf=rep(1,as.integer(max(group))),
				dfmax=as.integer(max(group))+1,pmax=min(dfmax*1.2,as.integer(max(group))),
				eps=1e-8, maxit=100,delta)
{
	#################################################################################	
	#	Design matrix setup, error checking
	this.call=match.call()
	loss = match.arg(loss)
    
    if(!is.matrix(x))
    stop("x has to be a matrix")
    
    if(any(is.na(x)))
    stop("Missing values in x not allowed!")
	
	y = drop(y)
	np=dim(x)
  	nobs=as.integer(np[1])
  	nvars=as.integer(np[2])
  	vnames=colnames(x)

  	if(is.null(vnames)) vnames=paste("V",seq(nvars),sep="")

	if(length(y)!=nobs) stop("x and y have different number of rows")

	if(!is.numeric(y)) 
	stop("The response y must be numeric. Factors must be converted to numeric")
	
    c1 = loss %in% c("logit","sqsvm","hsvm")
    c2 = any(y %in% c(-1,1) == FALSE)
    if(c1 && c2) 
    stop("Classification method requires the response y to be in {-1,1}")
	#################################################################################	
	#    group setup
	if(is.null(group)){group=1:nvars
	}else if (length(group)!=nvars) 
	stop("group length does not match the number of predictors in x")

    bn=as.integer(max(group))
    bs=as.integer(as.numeric(table(group)))

    if(!identical(as.integer(sort(unique(group))),as.integer(1:bn))) 
    stop("Groups must be consecutively numbered 1,2,3,...")

	ix=rep(NA,bn)
	iy=rep(NA,bn)
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
	#  get upper bound
	gamma <- rep(NA,bn)
	for(g in 1:bn) gamma[g] <- max(eigen(crossprod(x[,ix[g]:iy[g]]))$values)
	#################################################################################	
	#parameter setup
	if(missing(delta)) delta = 1
	if(delta<0) stop("delta must be non-negtive")
	delta=as.double(delta)
	if(missing(weights)) weights=rep(1,nobs)
    if(length(weights) != length(y))
    stop("length(weights) not equal length(y)")  
    if(any(weights < 0))
    stop("Negative weights not allowed")
    w = as.double(weights)
	if(length(pf)!=bn) 
	stop("The size of penalty factor must be same as the number of groups")
	maxit=as.integer(maxit)
	pf=as.double(pf)
	eps=as.double(eps)
	dfmax=as.integer(dfmax)
	pmax=as.integer(pmax)
	#################################################################################	
	#lambda setup
    nlam = as.integer(nlambda)
    if (is.null(lambda)) {
        if (lambda.factor >= 1) 
            stop("lambda.factor should be less than 1")
        flmin = as.double(lambda.factor)
        ulam = double(1)
    }
    else {
        #flmin=1 if user define lambda
        flmin = as.double(1)
        if (any(lambda < 0)) 
            stop("lambdas should be non-negative")
        ulam = as.double(rev(sort(lambda)))
        nlam = as.integer(length(lambda))
    }
	#################################################################################	
	# call Fortran core
	 fit=switch(loss,
		ls = ls(bn,bs,ix,iy,gamma,nobs,nvars,x,y,w,pf,dfmax,pmax,nlam,flmin,ulam,eps,maxit,vnames),
		hreg = hreg(delta,bn,bs,ix,iy,gamma,nobs,nvars,x,y,w,pf,dfmax,pmax,nlam,flmin,ulam,eps,maxit,vnames),
		logit = logit(bn,bs,ix,iy,gamma,nobs,nvars,x,y,w,pf,dfmax,pmax,nlam,flmin,ulam,eps,maxit,vnames),
		sqsvm = sqsvm(bn,bs,ix,iy,gamma,nobs,nvars,x,y,w,pf,dfmax,pmax,nlam,flmin,ulam,eps,maxit,vnames),
		hsvm = hsvm(delta,bn,bs,ix,iy,gamma,nobs,nvars,x,y,w,pf,dfmax,pmax,nlam,flmin,ulam,eps,maxit,vnames))
	#################################################################################	
	# output	
    if(is.null(lambda)) fit$lambda = lamfix(fit$lambda)
    fit$call = this.call
    class(fit) = c("gglasso", class(fit))
    fit
}