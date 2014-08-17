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
Yang, Y. and Zou, H. (2012), ``A Fast Unified Algorithm for Computing Group-Lasso Penalized Learning Problems,'' \emph{Statistics and Computing}. Accepted.\cr
BugReport: \url{http://code.google.com/p/gglasso/}\cr

Friedman, J., Hastie, T., and Tibshirani, R. (2010), ``Regularization paths for generalized
linear models via coordinate descent,'' \emph{Journal of Statistical Software}, 33, 1.\cr
\url{http://www.jstatsoft.org/v33/i01/}
}
\seealso{\code{\link{cv.gglasso}}.}
\examples{
# load gglasso library
library(gglasso)

# load data set
data(colon)

# define group index
group <- rep(1:20,each=5)

# 5-fold cross validation using group lasso 
# penalized logisitic regression
cv <- cv.gglasso(x=colon$x, y=colon$y, group=group, loss="logit",
pred.loss="misclass", lambda.factor=0.05, nfolds=5)

# make a CV plot
plot(cv)
}
\keyword{models}
\keyword{regression}
