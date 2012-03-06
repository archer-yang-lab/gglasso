source("bmd.exp.r")
source("err.r")
source("lamfix.r")
source("plot.bmd.r")
source("plot.bmd.exp.r")
dyn.load("bmdexp.so")

dl <- function(r,w)
{
	d=rep(0,length(r))
	for (i in 1:length(r))
	{
		if(r[i]>0) d[i]=2*w*r[i]
		else d[i]=2*r[i]
	}
	return (d)
}


set.seed(11)
x=matrix(rnorm(10*200),10,200) 
set.seed(11)
y=sample(c(-1,1),10,replace=T)
group<-rep(1:200,each=1)
nobs=nrow(x)
nvars=ncol(x)
#m0 <-bmd.exp(y=y,x=x,group=group,eps=1e-6)

bn=as.integer(max(group))
bs=as.integer(as.numeric(table(group)))
pf=rep(1,bn)

w=10
m1 <-bmd.exp(y=y,x=x,group=group,w=w,eps=1e-13,standardize=T,pf=pf)


one=rep(1, nobs)
meanx=drop(one %*% x)/nobs
x=scale(x, meanx, FALSE)

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
for (g in 1:bn)
{
	ind=ix[g]:iy[g]
	decomp <- qr(x[,ind])
	if(decomp$rank < bs[g]) stop("Block belonging to columns ",  ## Warn if block has not full rank
	paste(ind, collapse = ", ")," has not full rank! \n")
	x[,ind] <- qr.Q(decomp)
}

pf=pf*bn/sum(pf) 
B <- as.matrix(m1$beta)
MM <- matrix(0,bn,length(m1$lambda))
for (l in 1:length(m1$lambda))
{
	for (g in 1:bn)
	{	
		ind=(group==g)
		ri <- y-(x%*%B[,l]+m1$b0[l])
 		L = dl(ri,w)
		yxl <- t(x[,ind])%*%L
		yxlnorm <- sqrt(crossprod(yxl,yxl))
		Bnorm<-sqrt(crossprod(B[ind,l],B[ind,l]))
		
		if(sum(B[ind,l])!=0)
		{
			AA<- -yxl+  B[ind,l]*m1$lambda[l]*pf[g]*sqrt(bs[g])/Bnorm
			if(abs(sum(AA))<1e-10) MM[g,l] <- "."
			else{
				MM[g,l] <- "F"
				print(AA)
			} 
			
		}
		else
		{
			BB <- yxlnorm - pf[g] * m1$lambda[l] * sqrt(bs[g])
			if (BB<=0) MM[g,l] <- "."
			else {
					MM[g,l] <- "f"
					print(paste("this is",BB))
				}
		}
	}
}
print(MM)


