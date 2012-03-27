source("./source/gglasso.r")
source("./source/model.r")
source("./source/utilities.r")
source("./source/auxiliary.r")
dyn.load("./source/gglasso.so")


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


thr = 1e-6
loss = class(m1)[[2]]

KKT(b0 = m1$b0, m1$beta, y, x, m1$lambda, pf, group, thr, delta, loss = loss)