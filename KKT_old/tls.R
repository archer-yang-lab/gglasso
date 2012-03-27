source("./source/gglasso.r")
source("./source/model.r")
source("./source/utilities.r")
dyn.load("./source/gglasso.so")


dl <- function(r)
{
	d=rep(0,length(r))
	for (i in 1:length(r))
	{
		d[i] = r[i]
	}
	return (d)
}


set.seed(11)
n = 100
p = 200
x=matrix(rnorm(n*p),n,p) 
set.seed(11)
y=sample(c(-1,1),n,replace=T)
group<-rep(1:(p/5),each=5)
nobs=nrow(x)
nvars=ncol(x)
#m0 <-gglasso.ls(y=y,x=x,group=group,eps=1e-6)

bn=as.integer(max(group))
bs=as.integer(as.numeric(table(group)))
pf=rep(1,bn)
system.time(m1 <-gglasso(loss="ls",y=y,x=x,group=group,eps=1e-12,pf=pf))



B <- as.matrix(m1$beta)
for (l in 1:length(m1$lambda))
{
	for (g in 1:bn)
	{	
		ind=(group==g)
		ri <- y-(x%*%B[,l]+m1$b0[l])
 		L = dl(ri)
		yxl <- t(x[,ind])%*%L/nobs
		yxlnorm <- sqrt(crossprod(yxl,yxl))
		Bnorm<-sqrt(crossprod(B[ind,l],B[ind,l]))
		
		if(Bnorm!=0)
		{
			AA<- -yxl+  B[ind,l]*m1$lambda[l]*pf[g]/Bnorm
			if(abs(sum(AA)) >= 1e-5) print(abs(sum(AA)))
		}
		else
		{
			BB <- yxlnorm - pf[g] * m1$lambda[l]
			if (BB > 0) print(paste("this is",BB))
		}
	}
}


m1$df
