\name{cv.gglasso}
\alias{cv.gglasso}
\alias{cv.hsvmpath}
\alias{cv.sqsvmpath}
\alias{cv.logitpath}
\title{Cross-validation for gglasso}
\description{Does k-fold cross-validation for gglasso, produces a plot,
and returns a value for \code{lambda}. This function is modified based on the \code{cv} function from the \code{glmnet} package.}
\usage{
cv.gglasso(x, y, lambda, pred.loss, nfolds, foldid, delta,...)
}
\arguments{
		\item{x}{\code{x} matrix as in \code{\link{gglasso}}.}
		\item{y}{class label \code{y} as in \code{\link{gglasso}}.}
		\item{lambda}{optional user-supplied lambda sequence; default is
		\code{NULL}, and \code{\link{gglasso}} chooses its own sequence.}
		\item{nfolds}{number of folds - default is 5. Although \code{nfolds}
		can be as large as the sample size (leave-one-out CV), it is not
		recommended for large datasets. Smallest value allowable is \code{nfolds=3}.}
		\item{foldid}{an optional vector of values between 1 and \code{nfold}
		identifying what fold each observation is in. If supplied,
		\code{nfold} can be missing.}
		\item{pred.loss}{loss to use for cross-validation error. Valid options are:
		\itemize{
		\item \code{"loss"} Margin based loss function. 
		\item \code{"misclass"} gives misclassification error for classification.}
		Default is \code{"loss"}.}
		\item{delta}{parameter \eqn{\delta}{delta} only used in HHSVM for computing margin based loss function, only available with \code{pred.loss = "misclass"}.}
		\item{\dots}{other arguments that can be passed to gglasso.}

}
\details{The function runs \code{\link{gglasso}} \code{nfolds}+1 times; the
first to get the \code{lambda} sequence, and then the remainder to
compute the fit with each of the folds omitted. The average error and standard deviation over the
folds are computed. 
}
\value{an object of class \code{\link{cv.gglasso}} is returned, which is a
list with the ingredients of the cross-validation fit.
		\item{lambda}{the values of \code{lambda} used in the fits.}
		\item{cvm}{the mean cross-validated error - a vector of length
		\code{length(lambda)}.}
		\item{cvsd}{estimate of standard error of \code{cvm}.}
		\item{cvupper}{upper curve = \code{cvm+cvsd}.}
		\item{cvlower}{lower curve = \code{cvm-cvsd}.}
		\item{nzero}{number of non-zero coefficients at each \code{lambda}.}
		\item{name}{a text string indicating type of measure (for plotting
		purposes).}
		\item{gglasso.fit}{a fitted \code{\link{gglasso}} object for the full data.}
		\item{lambda.min}{The optimal value of \code{lambda} that gives minimum
		cross validation error \code{cvm}.}
		\item{lambda.1se}{The largest value of \code{lambda} such that error is
		within 1 standard error of the minimum.}
}

\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), "A Fast Unified Algorithm for Computing Group-Lasso Penalized Learning Problems," \emph{Journal of Computational and Graphical Statistics}. Under review.\cr
BugReport: \url{http://code.google.com/p/gglasso/}\cr
}

\seealso{\code{\link{gglasso}}, \code{\link{plot.cv.gglasso}}, \code{\link{predict.cv.gglasso}}, and \code{\link{coef.cv.gglasso}} methods.}
\examples{
data(FHT)
set.seed(2011)
cv=cv.gglasso(FHT$x,FHT$y,
lambda2=0.1,pred.loss="loss", nfolds=5)
plot(cv)
}
\keyword{models}
\keyword{regression}

