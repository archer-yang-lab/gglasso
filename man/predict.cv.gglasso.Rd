\name{predict.cv.gglasso}
\alias{predict.cv.gglasso}
\title{make predictions from a "cv.gglasso" object.}
\description{
This function makes predictions from a cross-validated gglasso model,
using the stored \code{"gglasso.fit"} object, and the optimal value
chosen for \code{lambda}.
}
\usage{
\method{predict}{cv.gglasso}(object, newx, s=c("lambda.1se","lambda.min"),...)
}
\arguments{
		\item{object}{fitted \code{\link{cv.gglasso}} object.}
		\item{newx}{matrix of new values for \code{x} at which predictions are
		to be made. Must be a matrix. See documentation for \code{predict.gglasso}.}
		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are required. Default is the value \code{s="lambda.1se"} stored
		on the CV object. Alternatively \code{s="lambda.min"} can be
		used. If \code{s} is numeric, it is taken as the value(s) of
		\code{lambda} to be used.}
		\item{\dots}{not used. Other arguments to predict. } }
\details{This function makes it easier to use the results of
  cross-validation to make a prediction.}
\value{The object returned depends the \dots argument which is passed on
to the \code{\link{predict}} method for \code{\link{gglasso}} objects.}

\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), "A Fast Unified Algorithm for Computing Group-Lasso Penalized Learning Problems," \emph{Journal of Computational and Graphical Statistics}. Under review.\cr
BugReport: \url{http://code.google.com/p/gglasso/}\cr
}

\seealso{\code{\link{cv.gglasso}}, and \code{\link{coef.cv.gglasso}} methods.}
\examples{
data(FHT)
set.seed(2011)
cv=cv.gglasso(FHT$x, FHT$y,
lambda2 = 1, pred.loss="misclass", 
lambda.factor=0.05,nfolds=5)
pre = predict(cv$gglasso.fit, newx = FHT$x, 
s = cv$lambda.1se, type = "class")
}
\keyword{models}
\keyword{regression}
