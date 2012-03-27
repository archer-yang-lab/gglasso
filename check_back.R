library(R.matlab)
source("source.R") # R code for generating the FHT data
require(MASS)   

load("env.rda")
source("./source/gglasso.r")
source("./source/plot.gglasso.r")
source("./source/loglik.r")
library(grplasso)
source("./source/model.r")
source("./source/utilities.r")
dyn.load("./source/gglasso.so")


system.time(m1 <-gglasso(loss="logit",y=yy,x=z,group=group,eps=1e-8,pf=sqrt(bs)))

path <- "./"
tmp <- readMat(file.path(path, "matlab.mat"))
beta = tmp$X
m=m1
m$beta = beta
m$group = group
m$lambda = m1$lambda

par(mfrow=c(2,1))
plot(m1)
plot(m)


pf=pf*bn/sum(pf) 
B <- as.matrix(m1$beta)
for (l in 1:length(m1$lambda))
{
	for (g in 1:bn)
	{	
		ind=(group==g)
		ri <- y*(x%*%B[,l]+m1$b0[l])
 		L = -1/(1+exp(ri))
		yxl <- t(x[,ind])%*%(L*y)/nobs
		yxlnorm <- sqrt(crossprod(yxl,yxl))
		Bnorm<-sqrt(crossprod(B[ind,l],B[ind,l]))
		
		if(Bnorm!=0)
		{
			AA<- yxl+  B[ind,l]*m1$lambda[l]*pf[g]/Bnorm
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