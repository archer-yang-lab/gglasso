source("source.R") # R code for generating the FHT data
require(MASS)   


source("./source/gglasso.r")
source("./source/plot.gglasso.r")
source("./source/loglik.r")
library(grplasso)
source("./source/model.r")
source("./source/utilities.r")
dyn.load("./source/gglasso.so")

n=100 # the number of observations
q=1000 # the number of predictors
rho = 0.5
x=genx2(n,q,rho)
z=genz(x,n,q)
y=genjerry(x,3)
y=as.vector(1*(runif(n)< 1/(1+exp(-y))))
yy = y
yy[y==0]=-1

write.table(yy,file="yy.txt",row.names=FALSE)
write.table(z,file="z.txt",row.names=FALSE,col.names=FALSE)


group<-rep(1:1000,each=3)
nobs=nrow(z)
nvars=ncol(z)

bn=as.integer(max(group))
bs=as.integer(as.numeric(table(group)))
pf=rep(1,bn)



system.time(m1 <-gglasso(loss="logit",y=yy,x=z,group=group,eps=1e-8,pf=sqrt(bs)))
lambda
write.table(lambda,file="lambda.txt",row.names=FALSE,col.names=FALSE)


zz = cbind(1,z)
ggroup = c(NA,group)
system.time(m2 <-grplasso(y=y,x=zz,index=ggroup, center = F, standardize = F, lambda = m1$lambda*n,control = grpl.control(tol=1e-8)))

m=m1
m$beta = m2$coef[-1,]
m$group = m2$index[-1]
m$lambda = m2$lambda

max(m1$beta - m2$coef[-1,])

    # m1_b0 <- m1$b0 - apply(meanx*m1$beta/normx,2,sum)
    # m1_beta <- m1$beta/normx

# l1 = loglik(x, yy, m1_beta, m1_b0)

par(mfrow=c(2,1))
plot(m1)
plot(m)
