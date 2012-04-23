\name{plot.cv.gcdnet}
\alias{plot.cv.gcdnet}
\title{plot the cross-validation curve produced by cv.gcdnet}
\description{
Plots the cross-validation curve, and upper and lower standard deviation
curves, as a function of the \code{lambda} values used. This function is modified based on the \code{plot.cv} function from the \code{glmnet} package.}
\usage{
\method{plot}{cv.gcdnet}(x, sign.lambda, ...)
}
\arguments{
	\item{x}{fitted \code{\link{cv.gcdnet}} object}
	\item{sign.lambda}{either plot against \code{log(lambda)} (default) or
	its negative if \code{sign.lambda=-1}.}
	\item{\dots}{other graphical parameters to plot}
}
\details{A plot is produced.}

\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), "An Efficient Algorithm for Computing The HHSVM and Its Generalizations," \emph{Journal of Computational and Graphical Statistics}. Accepted\cr
BugReport: \url{http://code.google.com/p/gcdnet/}\cr


Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}

\seealso{\code{\link{cv.gcdnet}}.}
\examples{
data(FHT)
set.seed(2011)
cv=cv.gcdnet(FHT$x, FHT$y,
lambda2 = 1, pred.loss="loss", nfolds=5)
plot(cv)
}
\keyword{models}
\keyword{regression}
