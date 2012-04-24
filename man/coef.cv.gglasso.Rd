\name{coef.cv.gglasso}
\alias{coef.cv.gglasso}
\title{get coefficients or make coefficient predictions from a "cv.gglasso" object.}
\description{
This function gets coefficients or makes coefficient predictions from a cross-validated gglasso model,
using the stored \code{"gglasso.fit"} object, and the optimal value
chosen for \code{lambda}.
}
\usage{
\method{coef}{cv.gglasso}(object,s=c("lambda.1se","lambda.min"),...)
}
\arguments{
	\item{object}{fitted \code{\link{cv.gglasso}} object.}
	\item{s}{value(s) of the penalty parameter \code{lambda} at which
	predictions are required. Default is the value \code{s="lambda.1se"} stored
	on the CV \code{object}, it is the largest value of \code{lambda} such that error is
	within 1 standard error of the minimum. Alternatively \code{s="lambda.min"} can be
	used, it is the optimal value of \code{lambda} that gives minimum
	cross validation error \code{cvm}. If \code{s} is numeric, it is taken as the value(s) of
	\code{lambda} to be used.}
	\item{\dots}{not used. Other arguments to predict. } 
}
\details{This function makes it easier to use the results of
cross-validation to get coefficients or make coefficient predictions.}
\value{The object returned depends the \dots argument which is passed on
to the \code{\link{predict}} method for \code{\link{gglasso}} objects.}
\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), "An Efficient Algorithm for Computing The HHSVM and Its Generalizations," \emph{Journal of Computational and Graphical Statistics}. Accepted\cr
BugReport: \url{http://code.google.com/p/gglasso/}\cr

Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}

\seealso{\code{\link{cv.gglasso}}, and \code{\link{predict.cv.gglasso}} methods.}
\examples{
data(FHT)
set.seed(2011)
cv=cv.gglasso(FHT$x, FHT$y, 
	lambda2 = 1,nfolds=5)
coef(cv,s="lambda.min")
}
\keyword{models}
\keyword{regression}