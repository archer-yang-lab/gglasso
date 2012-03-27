source("./source/gglasso.r")
source("./source/model.r")
source("./source/utilities.r")
dyn.load("./source/gglasso.so")

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
x=matrix(rnorm(50*200),50,200) 
set.seed(11)
y=sample(c(-1,1),50,replace=T)
group<-rep(1:40,each=5)
nobs=nrow(x)
nvars=ncol(x)


bn=as.integer(max(group))
bs=as.integer(as.numeric(table(group)))
delta=0.4
pf=rep(1,bn)
system.time(m1 <- gglasso(loss="hreg",y=y,x=x,group=group,eps=1e-12,pf=pf,delta=delta))


B <- as.matrix(m1$beta)
for (l in 1:length(m1$lambda))
{
	ri <- y-(x%*%B[,l]+m1$b0[l])
	for (g in 1:bn)
	{	
		ind=(group==g)
 		L = dl(ri,delta)
		yxl <- t(x[,ind])%*%L/nobs
		print(yxl)
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

