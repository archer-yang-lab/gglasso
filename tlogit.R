source("./source/gglasso.r")
source("./source/model.r")
source("./source/utilities.r")
source("./source/auxiliary.r")
dyn.load("./source/gglasso.so")


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
system.time(m1 <-gglasso(loss="logit",y=y,x=x,group=group,eps=1e-12,pf=pf))


thr = 1e-5
loss = class(m1)[[2]]

KKT(b0 = m1$b0, m1$beta, y, x, m1$lambda, pf, group, thr, delta, loss = loss)