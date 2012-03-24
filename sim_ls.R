source("gglasso.r")
source("plot.r")
source("loglik.r")
library(grpreg)
library(grplasso)
source("model.r")
source("utilities.r")
dyn.load("gglasso.so")
set.seed(1)
n = 1000
p = 20
x=matrix(rnorm(n*p),n,p) 
set.seed(1)
y=rnorm(n)



meanx <- apply(x,2,mean)
normx <- sqrt(apply((t(x)-meanx)^2,1,sum))/sqrt(n)
xx = scale(x,meanx,normx)


group<-rep(1:4,each=5)
nobs=nrow(x)
nvars=ncol(x)
#m0 <-gglasso.logit(y=y,x=x,group=group,eps=1e-6)

bn=as.integer(max(group))
bs=as.integer(as.numeric(table(group)))

#pf<-1:10
pf=rep(1,bn)
system.time(m1 <-gglasso(loss="ls",y=y,x=xx,group=group,eps=1e-8,pf=sqrt(bs)))
system.time(m2 <-grpreg(y=y,X=x,group=group,family="gaussian",penalty="gLasso",eps=1e-4,max.iter=1e8))



    m1_b0 <- m1$b0 - apply(meanx*m1$beta/normx,2,sum)
    m1_beta <- m1$beta/normx


l1 = loglik(x, yy, m1_beta, m1_b0)

par(mfrow=c(2,1))
plot(m1)
plot(m2)




pf=pf*bn/sum(pf) 
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
