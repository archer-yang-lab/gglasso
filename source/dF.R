dhreg <- function(r,delta)
{
	dl = 2.0 * sign(min(abs(r) , delta) , r)
	return (dl)
}


dhsvm <- function(r,delta)
{
	d=rep(0,length(r))
	for (i in 1:length(r))
	{
		if (r[i]>1) d[i]=0
		else 
		{
			if (r[i]<=(1-delta)) d[i]=-1
			else d[i]=(r[i]-1)/delta
		}
	}
	return (d)
}


dl <- function(r)
{
	d=rep(0,length(r))
	for (i in 1:length(r))
	{
		d[i] = r[i]
	}
	return (d)
}