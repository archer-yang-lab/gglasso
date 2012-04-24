\name{plot.cv.gglasso}
\alias{plot.cv.gglasso}
\title{plot the cross-validation curve produced by cv.gglasso}
\description{
Plots the cross-validation curve, and upper and lower standard deviation
curves, as a function of the \code{lambda} values used. This function is modified based on the \code{plot.cv} function from the \code{glmnet} package.}
\usage{
\method{plot}{cv.gglasso}(x, sign.lambda, ...)
}
\arguments{
	\item{x}{fitted \code{\link{cv.gglasso}} object}
	\item{sign.lambda}{either plot against \code{log(lambda)} (default) or
	its negative if \code{sign.lambda=-1}.}
	\item{\dots}{other graphical parameters to plot}
}
\details{A plot is produced.}

\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), "A Fast Unified Algorithm for Computing Group-Lasso Penalized Learning Problems," \emph{Journal of Computational and Graphical Statistics}. Under review.\cr
BugReport: \url{http://code.google.com/p/gglasso/}\cr
}
\seealso{\code{\link{cv.gglasso}}.}
\examples{
data(FHT)
set.seed(2011)
cv=cv.gglasso(FHT$x, FHT$y,
lambda2 = 1, pred.loss="loss", nfolds=5)
plot(cv)
}
\keyword{models}
\keyword{regression}
