nonzeroROW=function(beta)
{
	which=rep(FALSE,nrow(beta))
	for(i in 1:nrow(beta))
	{
		if(any(beta[i,]!=0)) which[i]=TRUE
	}
	res=(1:nrow(beta))[which]
	return(res)
}
