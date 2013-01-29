\name{coef.gglasso}
\alias{coef.gglasso}
\title{get coefficients or make coefficient predictions from an "gglasso" object.}
\description{
Computes the coefficients at the requested values for \code{lambda} from a fitted \code{\link{gglasso}} object.
}
\usage{
\method{coef}{gglasso}(object, s = NULL, ...)
}
\arguments{
		\item{object}{fitted \code{\link{gglasso}} model object.}
		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are required. Default is the entire sequence used to
		create the model.}
		\item{\dots}{not used. Other arguments to predict. } 
}

\details{
\code{s} is the new vector at which predictions are requested. If \code{s} is not in the lambda sequence used for fitting the model, the \code{coef} function will use linear interpolation to make predictions. The new values are interpolated using a fraction of coefficients from both left and right \code{lambda} indices. 
}
\value{The coefficients at the requested values for \code{lambda}.}


\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), ``A Fast Unified Algorithm for Computing Group-Lasso Penalized Learning Problems,'' submitted to \emph{Journal of Computational and Graphical Statistics}. Under review.\cr
BugReport: \url{http://code.google.com/p/gglasso/}\cr
}
\seealso{\code{\link{predict.gglasso}} method}
\examples{
# load gglasso library
library(gglasso)

# load data set
data(colon)

# define group index
group <- rep(1:20,each=5)

# fit group lasso
m1 <- gglasso(x=colon$x,y=colon$y,group=group,loss="logit")

# the coefficients at lambda = 0.01 and 0.02
coef(m1,s=c(0.01,0.02))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
