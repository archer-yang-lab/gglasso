\name{plot.gglasso}
\alias{plot.gglasso}
\title{Plot coefficients from a "gglasso" object}
\description{
Produces a coefficient profile plot of the coefficient paths for a
fitted \code{\link{gglasso}} object. This function is modified based on the \code{plot} function from the \code{glmnet} package.
}
\usage{
\method{plot}{gglasso}(x, xvar = c("norm", "lambda"), color = FALSE, label = FALSE, ...)
}
\arguments{
		\item{x}{fitted \code{\link{gglasso}} model}
		\item{xvar}{what is on the X-axis. \code{"norm"} plots against the
		L1-norm of the coefficients, \code{"lambda"} against the log-lambda
		sequence.}
		\item{color}{if \code{TRUE}, plot the curves with rainbow colors. \code{FALSE} is gray colors. Default is \code{FALSE}}
		\item{label}{if \code{TRUE}, label the curves with variable sequence numbers. Default is \code{FALSE}}
		\item{\dots}{other graphical parameters to plot}
}
\details{
A coefficient profile plot is produced.
}
\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), "A Fast Unified Algorithm for Computing Group-Lasso Penalized Learning Problems," \emph{Journal of Computational and Graphical Statistics}. Under review.\cr
BugReport: \url{http://code.google.com/p/gglasso/}\cr
}
\examples{
data(FHT)
m1 = gglasso(x=FHT$x,y=FHT$y)
par(mfrow=c(1,3))
plot(m1) # plots against the L1-norm of the coefficients
plot(m1,xvar="lambda",label=TRUE) # plots against the log-lambda sequence
plot(m1,color=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
