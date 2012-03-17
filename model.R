ls<-function(bn,bs,ix,iy,gamma,nobs,nvars,x,y,
				pf,dfmax,pmax,nlam,flmin,ulam,eps,maxit,vnames)
{
	#################################################################################	
	# call Fortran core
	fit=.Fortran("ls_f",bn,bs,ix,iy,gamma,
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
    outlist = getoutput(fit, maxit, pmax, nvars, vnames)
    outlist = c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
    class(outlist) = c("ls")
    outlist
}




logit<-function(bn,bs,ix,iy,gamma,nobs,nvars,x,y,
				pf,dfmax,pmax,nlam,flmin,ulam,eps,maxit,vnames)
{
	#################################################################################	
	# call Fortran core
	fit=.Fortran("log_f",bn,bs,ix,iy,gamma,
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
    outlist = getoutput(fit, maxit, pmax, nvars, vnames)
    outlist = c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
    class(outlist) = c("logit")
    outlist
}



hsvm<-function(delta,bn,bs,ix,iy,gamma,nobs,nvars,x,y,
				pf,dfmax,pmax,nlam,flmin,ulam,eps,maxit,vnames)
{
	#################################################################################	
	# call Fortran core
	fit=.Fortran("hsvm_f",delta,bn,bs,ix,iy,gamma,
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
    outlist = getoutput(fit, maxit, pmax, nvars, vnames)
    outlist = c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
    class(outlist) = c("hsvm")
    outlist
}


hreg<-function(delta,bn,bs,ix,iy,gamma,nobs,nvars,x,y,
				pf,dfmax,pmax,nlam,flmin,ulam,eps,maxit,vnames)
{
	#################################################################################	
	# call Fortran core
	fit=.Fortran("hreg_f",delta,bn,bs,ix,iy,gamma,
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
    outlist = getoutput(fit, maxit, pmax, nvars, vnames)
    outlist = c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
    class(outlist) = c("hreg")
    outlist
}