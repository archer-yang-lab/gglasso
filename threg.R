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
x=matrix(rnorm(100*200),100,200) 
set.seed(11)
y=sample(c(-1,1),100,replace=T)
group<-rep(1:40,each=5)
nobs=nrow(x)
nvars=ncol(x)


bn=as.integer(max(group))
bs=as.integer(as.numeric(table(group)))
delta=0.4
pf=rep(1,bn)
m1 <- gglasso(loss="hreg",y=y,x=x,group=group,eps=1e-8,pf=pf,delta=delta)


pf=pf*bn/sum(pf) 
B <- as.matrix(m1$beta)
for (l in 1:length(m1$lambda))
{
	for (g in 1:bn)
	{	
		ind=(group==g)
		ri <- y-(x%*%B[,l]+m1$b0[l])
 		L = dl(ri,delta)
		yxl <- t(x[,ind])%*%L/nobs
		yxlnorm <- sqrt(crossprod(yxl,yxl))
		Bnorm<-sqrt(crossprod(B[ind,l],B[ind,l]))
		if(Bnorm!=0)
		{
			AA<- -yxl+  B[ind,l]*m1$lambda[l]*pf[g]/Bnorm
			if(abs(sum(AA)) >= 1e-6) print(abs(sum(AA)))
		}
		else
		{
			BB <- yxlnorm - pf[g] * m1$lambda[l]
			if (BB > 0) print(paste("this is",BB))
		}
	}
}


