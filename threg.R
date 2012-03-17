source("gglasso.r")
source("model.r")
source("utilities.r")
dyn.load("gglasso.so")

dl <- function(r,delta)
{
	d=rep(0,length(r))
	for (i in 1:length(r))
	{
		if (r[i]< (-delta)) d[i]=-2*delta
		else if (r[i]>delta) d[i]=2*delta
		else d[i]=2*r[i]

	}
	return (d)
}


set.seed(11)
x=matrix(rnorm(10*200),10,200) 
set.seed(11)
y=sample(c(-1,1),10,replace=T)
group<-rep(1:40,each=5)
nobs=nrow(x)
nvars=ncol(x)


bn=as.integer(max(group))
bs=as.integer(as.numeric(table(group)))
delta=0.4
pf=rep(1,bn)
m1 <- gglasso(loss="hreg",y=y,x=x,group=group,eps=1e-13,standardize=T,pf=pf,delta=delta)


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
 		L = dl(ri,delta)
		yxl <- t(x[,ind])%*%L
		yxlnorm <- sqrt(crossprod(yxl,yxl))
		Bnorm<-sqrt(crossprod(B[ind,l],B[ind,l]))
		
		if(Bnorm!=0)
		{
			AA<- -yxl+  B[ind,l]*m1$lambda[l]*pf[g]*sqrt(bs[g])/Bnorm
			if(abs(sum(AA))<1e-9) MM[g,l] <- "."
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


