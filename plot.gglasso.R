plot.gglasso <- function(x, alpha=1, legend.loc, log.l=FALSE, ...)
  {
	xb = x$beta
	if(nrow(xb)==1){
		if (any(abs(xb) > 0)) {nonzeros = 1
		}else nonzeros = NULL
	}
	else{
		nonzeros <- which(apply(abs(xb),1,sum) > 0)
	}
    beta <- xb[nonzeros,,drop=FALSE]
    g <- as.numeric(as.factor(x$group[nonzeros]))
    p <- nrow(beta)
    l <- x$lambda
    n.g <- max(g)

    if (log.l)
      {
        l <- log(l)
        xlab <- expression(log(lambda))
      }
    else xlab <- expression(lambda)

    plot.args <- list(x=l, y=1:length(l), ylim=range(beta), xlab=xlab, ylab=expression(hat(beta)), type="n", xlim=rev(range(l)))
    new.args <- list(...)
    if (length(new.args))
      {
        new.plot.args <- new.args[names(new.args) %in% c(names(par()),names(formals(plot.default)))]
        plot.args[names(new.plot.args)] <- new.plot.args
      }
    do.call("plot", plot.args)

    line.args <- list(col=hcl(h=seq(0,360,len=(n.g+1)),l=70,c=100,alpha=alpha)[1:n.g],lwd=1+1.2^(-p/20),lty=1)
    if (length(new.args)) line.args[names(new.args)] <- new.args
    line.args$x <- l
    line.args$y <- t(beta)
    line.args$col <- rep(line.args$col,table(g))
    do.call("matlines",line.args)

    abline(h=0,lwd=line.args$lwd)
    
    if(!missing(legend.loc))
      {
        legend.args <- list(col=hcl(h=seq(0,360,len=(n.g+1)),l=70,c=100,alpha=alpha)[1:n.g],lwd=1+1.2^(-p/20),lty=1,legend=unique(g))
        if (length(new.args))
          {
            new.legend.args <- new.args[names(new.args) %in% names(formals(legend))]
            legend.args[names(new.legend.args)] <- new.legend.args
          }
        legend.args$x <- legend.loc
        do.call("legend",legend.args)
      }
  }
