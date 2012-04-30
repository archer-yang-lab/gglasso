\name{plot.gglasso}
\alias{plot.gglasso}
\title{Plot solution paths from a "gglasso" object}
\description{
Produces a coefficient profile plot of the coefficient paths for a
fitted \code{\link{gglasso}} object.
}
\usage{
\method{plot}{gglasso}(x, group = FALSE, log.l = TRUE, ...)
}
\arguments{
		\item{x}{fitted \code{\link{gglasso}} model}
		\item{group}{what is on the Y-axis. Plot the norm of each group
		 if \code{TRUE}. Plot each coefficient if \code{FALSE}.}
		\item{log.l}{what is on the X-axis. Plot against the log-lambda
		sequence if \code{TRUE}. Plot against the lambda sequence if \code{FALSE}.}
		\item{\dots}{other graphical parameters to plot}
}
\details{
A coefficient profile plot is produced.
}
\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), ``A Fast Unified Algorithm for Computing Group-Lasso Penalized Learning Problems,'' submitted to \emph{Journal of Computational and Graphical Statistics}. Under review.\cr
BugReport: \url{http://code.google.com/p/gglasso/}\cr
}
\examples{
# load gglasso library
library(gglasso)

# load data set
data(bardet)

# define group index
group <- rep(1:200,each=5)

# fit group lasso
m1 <- gglasso(x=bardet$x,y=bardet$y,group=group,loss="ls")

# make plots
par(mfrow=c(1,3))
plot(m1) # plots the coefficients against the log-lambda sequence 
plot(m1,group=TRUE) # plots group norm against the log-lambda sequence 
plot(m1,log.l=FALSE) # plots against the lambda sequence
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
