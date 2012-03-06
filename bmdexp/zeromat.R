zeromat=function(nvars,nalam,vnames,stepnames){
  ca=rep(0,nalam)
  ia=seq(nalam+1)
  ja=rep(1,nalam)
  dd=c(nvars,nalam)
  new("dgCMatrix", Dim = dd,
      Dimnames = list(vnames,stepnames),
      x = as.vector(ca),
      p = as.integer(ia - 1), i = as.integer(ja - 1))
}
