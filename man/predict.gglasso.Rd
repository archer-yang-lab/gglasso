\name{predict.gglasso}
\alias{predict.gglasso}
\alias{predict.hsvmpath}
\alias{predict.sqsvmpath}
\alias{predict.logitpath}
\title{make predictions from a "gglasso" object.}

\description{
Similar to other predict methods, this functions predicts fitted values and class labels from a fitted \code{\link{gglasso}} object.
}

\usage{
\method{predict}{gglasso}(object, newx, s = NULL,
type=c("class","link"), ...)
}

\arguments{
		\item{object}{fitted \code{\link{gglasso}} model object.}
		\item{newx}{matrix of new values for \code{x} at which predictions are
		to be made. Must be a matrix.}
		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are required. Default is the entire sequence used to
		create the model. }
		\item{type}{
		type of prediction required. 
		\itemize{
				\item Type \code{"link"} gives the
				linear predictors for classification problems.
				\item Type \code{"class"} produces the class label corresponding to the maximum probability.}
		}
		\item{\dots}{Not used. Other arguments to predict. } }
\details{\code{s} is the new vector at which predictions are requested. If \code{s} is not in the lambda sequence used for fitting the model, the \code{predict} function will use linear interpolation to make predictions. The new values are interpolated using a fraction of predicted values from both left and right \code{lambda} indices. }
\value{The object returned depends on type.}


\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), "A Fast Unified Algorithm for Computing Group-Lasso Penalized Learning Problems," \emph{Journal of Computational and Graphical Statistics}. Under review.\cr
BugReport: \url{http://code.google.com/p/gglasso/}\cr
}
\seealso{\code{\link{coef}} method}
\examples{
data(FHT)
m1 = gglasso(x=FHT$x,y=FHT$y)
print(predict(m1,type="class",newx=FHT$x[2:5,]))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
