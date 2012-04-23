\name{gcdnet}
\alias{gcdnet}
\title{Fits the regularization paths for large margin classifiers}
\description{Fits a regularization path for large margin classifiers at a sequence of regularization parameters lambda.}
\usage{
gcdnet(x, y, nlambda = 100, 
		method = c("hhsvm", "logit", "sqsvm"),
		lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04), lambda = NULL, lambda2 = 0, 
		pf = rep(1, nvars), exclude, dfmax = nvars + 1, pmax = min(dfmax * 1.2, 
	    nvars), standardize = TRUE, eps = 1e-8, maxit = 1e6, delta = 2)
}
\arguments{
		\item{x}{matrix of predictors, of dimension \eqn{N \times p}{N*p}; each row is an observation vector.}

		\item{y}{response variable. This argument should be a two-level factor for classification.}
		\item{nlambda}{the number of \code{lambda} values - default is 100.}
		\item{method}{a character string specifying the loss function to use, valid options are:
		\itemize{
		\item \code{"hhsvm"} Huberized squared hinge loss,
		\item \code{"sqsvm"} Squared hinge loss,
		\item \code{"logit"} logistic loss.
		}Default is \code{"loss"}.}
		\item{lambda.factor}{The factor for getting the minimal lambda in \code{lambda} sequence, where \code{min(lambda) = lambda.factor * max(lambda)}.  \code{max(lambda)} is the smallest value of \code{lambda} for which all coefficients are zero. The default depends on the relationship between \eqn{N} (the number of rows in the matrix of predictors) and \eqn{p} (the number of predictors). If \eqn{N > p}, the default is \code{0.0001},
		close to zero.  If \eqn{N<p}, the default is \code{0.01}.
		A very small value of \code{lambda.factor} will lead to a saturated fit. It takes no effect if there is user-defined \code{lambda} sequence.} 
		\item{lambda}{a user supplied \code{lambda} sequence. Typically, by leaving this option unspecified users can have 
		the program compute its own \code{lambda} sequence based on
		\code{nlambda} and \code{lambda.factor}. Supplying a value of
		\code{lambda} overrides this. It is better to supply
		a decreasing sequence of \code{lambda} values than a single (small) value, if not, the program will sort user-defined \code{lambda} sequence in decreasing order automatically.}

		\item{lambda2}{regularization parameter \eqn{\lambda_2}{lambda2} for the quadratic penalty of the 
		coefficients.}

		\item{pf}{penalty factor of length \eqn{p}. Separate penalty weights can be applied to each coefficient of \eqn{\beta}{beta} to allow
		differential shrinkage. Can be 0 for some variables, which implies
		no shrinkage, and results in that variable always being included in the
		model. Default is 1 for all variables (and implicitly infinity for
		variables listed in \code{exclude}).}

		\item{exclude}{indices of variables to be excluded from the
		model. Default is none. Equivalent to an infinite penalty factor.}

		\item{dfmax}{limit the maximum number of variables in the
		model. Useful for very large \eqn{p}, if a partial path is desired. Default is \eqn{p+1}.}

		\item{pmax}{limit the maximum number of variables ever to be nonzero. For example once \eqn{\beta} enters the model, no matter how many times it exits or re-enters model through the path, it will be counted only once. Default is \code{min(dfmax*1.2,p)}.}

		\item{standardize}{logical flag for variable standardization, prior to
		fitting the model sequence. If \code{TRUE}, x matrix is normalized such that sum squares of each column \eqn{\sum^N_{i=1}x_{ij}^2/N=1}{<Xj,Xj>/N=1}. Note that x is always centered (i.e. \eqn{\sum^N_{i=1}x_{ij}=0}{sum(Xj)=0}) no matter \code{standardize} is \code{TRUE} or \code{FALSE}. The coefficients are always returned on
		the original scale. Default is is \code{TRUE}.}

		\item{eps}{convergence threshold for coordinate majorization descent. Each inner
		coordinate majorization descent loop continues until the relative change in any
		coefficient (i.e. \eqn{\max_j(\beta_j^{new}-\beta_j^{old})^2}{max(j)(beta_new[j]-beta_old[j])^2}) is less than \code{eps}. For HHSVM i.e. \code{method="hhsvm"}, it is \eqn{\frac{2}{\delta}\max_j(\beta_j^{new}-\beta_j^{old})^2}{2*max(j)(beta_new[j]-beta_old[j])^2/delta}. Defaults value is \code{1e-8}.}

		\item{maxit}{maximum number of outer-loop iterations allowed at fixed lambda value. Default is 1e6. If models do not converge, consider increasing \code{maxit}.}

		\item{delta}{the parameter \eqn{\delta}{delta} in HHSVM model. Default is 2.}
}

\details{
Note that the objective function in \code{gcdnet} is
\deqn{Loss(y, X, beta))/N + lambda1 * |beta| + 0.5 * lambda2 * beta^2}
where the penalty is a combination of L1 and L2 term. Users can specify the loss function to use, options include Huberized squared hinge loss, Squared hinge loss and logistic regression. Users can also tweak the penalty by choosing different \eqn{lambda2} and penalty factor. 

For computing speed reason, if models are not converging or running slow, consider increasing \code{eps}, decreasing
\code{nlambda}, or increasing \code{lambda.factor} before increasing
\code{maxit}.
}


\value{
An object with S3 class \code{\link{gcdnet}}.
		\item{call}{the call that produced this object}
		\item{b0}{intercept sequence of length \code{length(lambda)}}
		\item{beta}{a \code{p*length(lambda)} matrix of coefficients, stored as a sparse matrix (\code{dgCMatrix} class, the standard class for sparse numeric matrices in the \code{Matrix} package.). To convert it into normal type matrix use \code{as.matrix()}.}
		\item{lambda}{the actual sequence of \code{lambda} values used}
		\item{df}{the number of nonzero coefficients for each value of
		\code{lambda}.}
		\item{dim}{dimension of coefficient matrix (ices)}
		\item{npasses}{total number of iterations (the most inner loop) summed over all lambda values}
		\item{jerr}{error flag, for warnings and errors, 0 if no error.}
}

\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), "An Efficient Algorithm for Computing The HHSVM and Its Generalizations," \emph{Journal of Computational and Graphical Statistics}. Accepted\cr
BugReport: \url{http://code.google.com/p/gcdnet/}\cr
}


\seealso{\code{plot.gcdnet}}
\examples{
data(FHT)
# solution paths for the HHSVM
m1 <- gcdnet(x=FHT$x,y=FHT$y,delta=1,lambda2=1,method="hhsvm")
plot(m1)
# solution paths for the penalized SVM with the squared hinge loss
m2 <- gcdnet(x=FHT$x,y=FHT$y,lambda2=0.1,method="sqsvm")
plot(m2)
# solution paths for the penalized logistic regression
m3 <- gcdnet(x=FHT$x,y=FHT$y,lambda2=0.01,method="logit")
plot(m3)
}
\keyword{models}
\keyword{regression}
