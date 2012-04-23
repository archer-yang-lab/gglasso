\name{predict.cv.gcdnet}
\alias{predict.cv.gcdnet}
\title{make predictions from a "cv.gcdnet" object.}
\description{
This function makes predictions from a cross-validated gcdnet model,
using the stored \code{"gcdnet.fit"} object, and the optimal value
chosen for \code{lambda}.
}
\usage{
\method{predict}{cv.gcdnet}(object, newx, s=c("lambda.1se","lambda.min"),...)
}
\arguments{
		\item{object}{fitted \code{\link{cv.gcdnet}} object.}
		\item{newx}{matrix of new values for \code{x} at which predictions are
		to be made. Must be a matrix. See documentation for \code{predict.gcdnet}.}
		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are required. Default is the value \code{s="lambda.1se"} stored
		on the CV object. Alternatively \code{s="lambda.min"} can be
		used. If \code{s} is numeric, it is taken as the value(s) of
		\code{lambda} to be used.}
		\item{\dots}{not used. Other arguments to predict. } }
\details{This function makes it easier to use the results of
  cross-validation to make a prediction.}
\value{The object returned depends the \dots argument which is passed on
to the \code{\link{predict}} method for \code{\link{gcdnet}} objects.}

\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), "An Efficient Algorithm for Computing The HHSVM and Its Generalizations," \emph{Journal of Computational and Graphical Statistics}. Accepted\cr
BugReport: \url{http://code.google.com/p/gcdnet/}\cr


Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}


\seealso{\code{\link{cv.gcdnet}}, and \code{\link{coef.cv.gcdnet}} methods.}
\examples{
data(FHT)
set.seed(2011)
cv=cv.gcdnet(FHT$x, FHT$y,
lambda2 = 1, pred.loss="misclass", 
lambda.factor=0.05,nfolds=5)
pre = predict(cv$gcdnet.fit, newx = FHT$x, 
s = cv$lambda.1se, type = "class")
}
\keyword{models}
\keyword{regression}
