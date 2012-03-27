genx2=function(n,q,rho){
#    generate x's multivariate normal with equal corr rho
 if(abs(rho)<1){
  beta=sqrt(rho/(1-rho))
 x0=matrix(rnorm(n*q),ncol=q)
 z=rnorm(n)
 x=beta*matrix(z,nrow=n,ncol=q,byrow=F)+x0
 }
 if(abs(rho)==1){ x=matrix(z,nrow=n,ncol=q,byrow=F)}

return(x)
}

genz=function(x,n,q){
z = matrix(NA,n,3*q) 
for(k in 1:q){
	i1 = 3*k - 2
	i2 = 3*k - 1
	i3 = 3*k
	z[,i1] = 2*x[,k]/3
	z[,i2] = -x[,k]^2
	z[,i3] = x[,k]^3/3
}
return(z)
}

genjerry=
function(x,snr){
# generate data according to Friedman's setup
n=nrow(x)
q=ncol(x)
b=((-1)^(1:q))*exp(-(2*(1:q)-1)/20)
f=x%*%b
Ystar=2*f/3-f^2+f^3/3
e=rnorm(n)
k=sqrt(var(Ystar)/(snr*var(e)))
y=Ystar+k*e
return(y)
}

