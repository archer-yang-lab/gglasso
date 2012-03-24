source("gglasso.r")
source("plot.gglasso.r")
source("loglik.r")
library(grplasso)
source("model.r")
source("utilities.r")
dyn.load("gglasso.so")
set.seed(1)
n = 100
p = 2000
x=matrix(rnorm(n*p),n,p) 
set.seed(1)
y=sample(c(-1,1),n,replace=T)
yy = y
yy[y == -1] = 0


meanx <- apply(x,2,mean)
normx <- sqrt(apply((t(x)-meanx)^2,1,sum))/sqrt(n)
xx = scale(x,meanx,normx)


group<-rep(1:20,each=100)
nobs=nrow(x)
nvars=ncol(x)
#m0 <-gglasso.logit(y=y,x=x,group=group,eps=1e-6)

bn=as.integer(max(group))
bs=as.integer(as.numeric(table(group)))

#pf<-1:10
pf=rep(1,bn)
system.time(m1 <-gglasso(loss="logit",y=y,x=xx,group=group,eps=1e-8,pf=sqrt(bs)))

xxx = cbind(1,xx)
ggroup = c(NA,group)


# lambda = lambdamax(y=yy,x=xxx,index=ggroup, center = F, standardize = F)* 0.5^(0:100)
system.time(m3 <-grplasso(y=yy,x=xxx,index=ggroup, center = F, standardize = F, lambda = m1$lambda*n,control = grpl.control(tol=1e-8)))

# lambda[1]/m1$lambda[1]

max(m1$beta - m3$coef[-1,])


    # m1_b0 <- m1$b0 - apply(meanx*m1$beta/normx,2,sum)
    # m1_beta <- m1$beta/normx


# l1 = loglik(x, yy, m1_beta, m1_b0)

par(mfrow=c(2,1))
plot(m1)
plot(m3)
