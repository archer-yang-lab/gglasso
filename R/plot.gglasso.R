plot.gglasso <- function(x, alpha = 1, type = c("indiv", "group"), 
    legend.loc, log.l = FALSE, ...) {
    type <- match.arg(type)
    xb <- x$beta
    if (nrow(xb) == 1) {
        if (any(abs(xb) > 0)) {
            nonzeros <- 1
        } else nonzeros <- NULL
    } else {
        nonzeros <- which(apply(abs(xb), 1, sum) > 0)
    }
    tmp <- xb[nonzeros, , drop = FALSE]
    g <- as.numeric(as.factor(x$group[nonzeros]))
    p <- nrow(tmp)
    l <- x$lambda
    n.g <- max(g)
    
    switch(type, group = {
        bs <- as.integer(as.numeric(table(g)))
        ix <- rep(NA, n.g)
        iy <- rep(NA, n.g)
        j <- 1
        for (g in 1:n.g) {
            ix[g] <- j
            iy[g] <- j + bs[g] - 1
            j <- j + bs[g]
        }
        beta <- matrix(NA, n.g, length(l))
        for (g in 1:n.g) {
            crossp <- apply(tmp[ix[g]:iy[g], ], 2, crossprod)
            beta[g, ] <- sqrt(crossp)
        }
    }, indiv = {
        beta <- tmp
    })
    
    if (log.l) {
        l <- log(l)
        xlab <- "Log Lambda"
    } else xlab <- "Lambda"
    
    plot.args <- list(x = l, y = 1:length(l), ylim = range(beta), xlab = xlab, 
        ylab = "Coefficients", type = "n", xlim = range(l))
    new.args <- list(...)
    if (length(new.args)) {
        new.plot.args <- new.args[names(new.args) %in% c(names(par()), names(formals(plot.default)))]
        plot.args[names(new.plot.args)] <- new.plot.args
    }
    do.call("plot", plot.args)
    line.args <- list(col = gray.colors(n.g + 1, start = 0.05, end = 0.7, gamma = 2.2)[1:n.g], 
        lwd = 1 + 1.2^(-p/20), lty = 1)
    
    if (length(new.args)) 
        line.args[names(new.args)] <- new.args
    line.args$x <- l
    line.args$y <- t(beta)
    line.args$col <- rep(line.args$col, table(g))
    do.call("matlines", line.args)
    
    abline(h = 0, lwd = line.args$lwd)
    
    if (!missing(legend.loc)) {
        legend.args <- list(col = hcl(h = seq(0, 360, len = (n.g + 1)), l = 70, 
            c = 100, alpha = alpha)[1:n.g], lwd = 1 + 1.2^(-p/20), lty = 1, 
            legend = unique(g))
        if (length(new.args)) {
            new.legend.args <- new.args[names(new.args) %in% names(formals(legend))]
            legend.args[names(new.legend.args)] <- new.legend.args
        }
        legend.args$x <- legend.loc
        do.call("legend", legend.args)
    }
} 
