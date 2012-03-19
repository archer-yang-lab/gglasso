source("gglasso.r")
source("model.r")
source("utilities.r")
dyn.load("gglasso.so")
set.seed(1)
x=matrix(rnorm(100*200),100,200) 
set.seed(1)
y=sample(c(-1,1),100,replace=T)
group<-rep(1:40,each=5)
nobs=nrow(x)
nvars=ncol(x)
#m0 <-gglasso.logit(y=y,x=x,group=group,eps=1e-6)

bn=as.integer(max(group))
bs=as.integer(as.numeric(table(group)))

#pf<-1:10
pf=rep(1,bn)
system.time(m1 <-gglasso(loss="sqsvm",y=y,x=x,group=group,eps=1e-8,pf=pf))


pf=pf*bn/sum(pf) 
B <- as.matrix(m1$beta)
for (l in 1:length(m1$lambda))
{
	for (g in 1:bn)
	{	
		ind=(group==g)
		ri <- y*(x%*%B[,l]+m1$b0[l])
 		L= -2*ifelse(ri<=1,(1-ri),0)
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


