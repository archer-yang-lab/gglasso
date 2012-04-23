coef.gglasso = function(object, s = NULL, ...) {
    b0 = t(as.matrix(object$b0))
    rownames(b0) = "(Intercept)"
    nbeta = rbind2(b0, object$beta)
    if (!is.null(s)) {
        vnames = dimnames(nbeta)[[1]]
        dimnames(nbeta) = list(NULL, NULL)
        lambda = object$lambda
        lamlist = lambda.interp(lambda, s)
        nbeta = nbeta[, lamlist$left, drop = FALSE] * lamlist$frac + 
            nbeta[, lamlist$right, drop = FALSE] * (1 - lamlist$frac)
        dimnames(nbeta) = list(vnames, paste(seq(along = s)))
    }
    return(nbeta)
} 


predict.gglasso = function(object, newx, s = NULL, type = c("class", "link"), ...){
	type = match.arg(type)
	loss = class(object)[[2]]
    b0 = t(as.matrix(object$b0))
    rownames(b0) = "(Intercept)"
    nbeta = rbind2(b0, object$beta)
    if (!is.null(s)) {
        vnames = dimnames(nbeta)[[1]]
        dimnames(nbeta) = list(NULL, NULL)
        lambda = object$lambda
        lamlist = lambda.interp(lambda, s)
        nbeta = nbeta[, lamlist$left, drop = FALSE] * lamlist$frac + 
            nbeta[, lamlist$right, drop = FALSE] * (1 - lamlist$frac)
        dimnames(nbeta) = list(vnames, paste(seq(along = s)))
    }
    nfit = as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
	if(loss %in% c("logit","sqsvm","hsvm")){
		switch(type, link = nfit, class = ifelse(nfit > 0, 1, -1))
	}
	else{nfit}
}