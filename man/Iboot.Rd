\name{Iboot}
\alias{Iboot}
\title{Iterated bootstrap tests and confidence sets for a \eqn{p}-dimensional parameter}
\description{
   This function provides both simple and re-calibrated bootstrap critical values for tests and confidence sets for a \eqn{p}-dimensional parameter, obtained by a single weighted bootstrap iteration. The considered statistic is the unstudentised version of the Rao statistic as proposed in Lunardon (2013). Details about the implemented algorithm can be found in Nankervis (2005) and Lunardon (2013). 
}
\usage{
Iboot(gradient, B=500, M=500, kB=0.01, alpha=c(0.1, 0.05), control.optim=list(), seed)
}
\arguments{
   \item{gradient}{An object of class \code{data.matrix}, or coercible to that class, of dimension \code{n x p}, where \code{n} is the sample size and \code{p} is the dimension of the parameter, containing the contributions of an estimating function for the parameter evaluated at the desired value. See ``Details''.}
   \item{B}{An integer indicating the number of outer level bootstrap replications.}
   \item{M}{An integer indicating the number of inner level bootstrap replications. If equal to 0 then single bootstrap is performed.}
   \item{kB}{The proportion of convex hull condition failures in the outer level that the user allows before stopping the algorithm. See ``Details''.}
   \item{alpha}{A vector specifying the desired nominal level(s) of test and confidence set.}
	\item{control.optim}{A list of optional control parameters as passed to the optimisation routine \code{\link{optim}} for \code{method="L-BFGS-B"}. See ``Details''.} 
   \item{seed}{A single value, interpreted as an integer, recommended to specify seeds.}
}

\details{
\code{Iboot} performs a single weighted bootstrap iteration in order to provide re-calibrated critical values for tests and confidence sets based on the unstudentised version of the Rao statistic for a parameter \eqn{\theta\in R^p,\,p\geq 1}. 

Denoted with \eqn{y=(y_1,\dots,y_n)} an i.i.d. sample of size \eqn{n}, where \eqn{y_i\in R^d,\,d\geq 1}, and with \eqn{U(\theta)=\sum_{i=1}^n u(\theta;y_i)} an (unbiased) estimating function for \eqn{\theta}, the unstudentised version of the Rao statistic is defined as
\deqn{
W_{us}(\theta)=n^{-1} U(\theta)^\top U(\theta).
}

The function \code{Iboot} implements the algorithm outlined in Lunardon (2013) that differs from the Nankervis's one for two sided confidence intervals because bootstrap is carried out in a weighted fashion. In the outer level bootstrap versions of \eqn{W_{us}(\theta)}, \eqn{W^{b}_{us}(\theta),\,b=1,\dots,B}, are obtained by resampling according to the \eqn{n}-dimensional vector of weights \eqn{w(\theta)=(w_1(\theta),\dots,w_n(\theta))} associated to each contribution \eqn{u(\theta;y_i)}. The functional form of \eqn{w_i(\theta)} is provided by Owen's empirical likelihood formulation, that is
\deqn{
w_i(\theta)=n^{-1}(1+\lambda(\theta)^\top u(\theta;y_i))^{-1},
}
where \eqn{\lambda(\theta)\in R^p} is the Lagrange multiplier (see Owen, 1990). Instead, in the inner level, bootstrap versions of \eqn{W^{b}_{us}(\theta)} are obtained by resampling according to \eqn{w^b(\theta)}, computed as before by using the outer level contributions \eqn{u(\theta; y^b_i)}.

As weighted bootstrap might entail computational difficulties, i.e. the convex hull condition (see Owen, 1990) might not be satisfied in all \code{B} bootstrap replications, some precautions have been taken. In particular, the algorithm checks if the convex hull condition is satisfied for both the observed contributions \eqn{u(\theta;y_i)} and for each of the \code{B} bootstrap resamplings in the outer level, i.e. \eqn{u(\theta;y^{b}_i)}. 

The algorithm stops immediately if the convex hull condition fails for \eqn{u(\theta;y_i)} whereas concludes after reaching \code{B x kB} convex hull condition failures in the outer level. In particular, in the former situation the observed value of the statistic is returned only, whereas in the latter both the observed value and the bootstrap distribution (based on a smaller number of resamplings than \code{B}) of the statistic are supplied.

Numerical optimisation of the objective function that supplies \eqn{\lambda(\theta)} relies on the foreign function \code{lbfgsb} called from R by \code{optim} with \code{method="L-BFGS-B"}. However, it is not possible to set bounds on the variables by specifying \code{lower} and \code{upper}, instead the optional parameters that can be set in \code{control.optim} are

  \describe{
  		\item{\code{trace}}{Non-negative integer. If positive, tracing information on the progress of the optimisation is produced. Higher values
          may produce more tracing information (up to six levels of tracing). Defaults to 0.}
  		\item{\code{maxit}}{The maximum number of iterations. Defaults to 2e4.}
  		\item{\code{pgtol}}{It is a tolerance on the projected gradient in the current search direction. Defaults to \code{sqrt(.Machine$double.eps)}.}
  		\item{\code{REPORT}}{The frequency of reports  if \code{trace} is positive. Defaults to every 10 iterations.}
  		\item{\code{lmm}}{An integer giving the number of BFGS updates retained. Defaults to 5.}
  		\item{\code{factr}}{Controls the convergence. Convergence occurs when the reduction in the objective is within this factor of the machine tolerance. Default is 1e7, that is a tolerance of about 1e-8.}
  }
}

\value{
A list of class ``Iboot'' containing:
  \item{Call}{The matched call.}
  \item{oss}{ The observed value of the statistic.}
  \item{boot}{The bootstrap distribution of the statistic (sorted into descending order).}
  \item{map}{The re-calibrated \code{1-alpha} nominal levels.}
  \item{boot.quant}{The \code{1-alpha} level bootstrap quantile(s).}
  \item{recalib.quant}{The re-calibrated \code{1-alpha} level bootstrap quantile(s).}
  \item{fails.outer}{Actual proportion of failures in the outer level.}
  \item{failure}{An error code indicating wheter the algorithm has succeeded:
  \describe{
  		\item{\code{0}}{indicates that the algorithm have reached the end successfully (see \code{fails.outer} to check the actual proportion of failures).}
  		\item{\code{1}}{convex hull condition not satisfied by using the original contributions \eqn{u(\theta;y_i)}.}
  		\item{\code{2}}{the actual proportion of convex hull condition failures in the outer level has exceeded \code{kB}.}
  }}
}

\note{
S3 print and summary methods are associated to objects of class ``Iboot''.
}

\references{
Lunardon, N. (2013). Prepivoting composite score statistics by weighted bootstrap iteration. E-print: arXiv/1301.7026.

Nankervis, J. (2005). Computational algorithms for double bootstrap confidence intervals. \emph{Computational statistics & data analysis}, \bold{45}, 461--475.

Owen, A. (1990). Empirical likelihood ratio confidence regions. \emph{The Annals of Statistics}, \bold{18}, 90--120.


}

\seealso{
\code{\link[simpleboot:one.boot]{one.boot}}, \code{\link{boot}}, \code{\link{optim}} 
}
\examples{

####Example 1: mean of a normal with known scale
n <- 20
mu <- 1

set.seed(1)
##contributions obtained from the score function
gr <- rnorm(n, mu) - mu

OBJ.Ib <- Iboot(gradient=gr, B=500, M=500, kB=0.01, alpha=c(0.1, 0.05, 0.01), seed=1)

##critical values for testing H0: mu=1, H1: mu!=1
OBJ.Ib
summary(OBJ.Ib)

#####exceeded number of convex hull condition failures in the outer level (kB=0)
OBJ.Ib <- Iboot(gradient=gr, B=500, M=500, kB=0, alpha=c(0.1, 0.05, 0.01), seed=1)
OBJ.Ib
summary(OBJ.Ib)

\dontrun{
####Example 2: example 1 of Lunardon (2013)
n <- 20
q <- 10

##parameter
mu <- 0
sig2 <- 1
rho <- 0.5
theta <- c(mu, sig2, rho)

##function to compute the pairwise score contributions
pscore_theta <- function(theta, data) 
{
	n <- nrow(data)
	q <- ncol(data)
	mu <- theta[1]
	sig2 <- theta[2]
	rho <- theta[3]
	A <- matrix(rho, q, q)
	diag(A) <- -(q-1)
	A <- A/((1-rho^2)*sig2^2)
	B <- matrix(-(1+rho^2), q, q)
	diag(B) <- 2*rho*(q-1)
	B <- B/(sig2*(1-rho^2)^2)
	x_dot <- apply(data, 1, sum)
	p_mu <- ((q-1)/(sig2*(1+rho)))*(x_dot-q*mu)
	p_sig2 <- -0.5*apply(((data-mu)%*%A)*(data-mu),1,sum)-q*(q-1)*0.5/sig2
	p_rho <- q*(q-1)*rho*0.5/(1-rho^2)-0.5*apply(((data-mu)%*%B)*(data-mu),1,sum)
	RES <- cbind(p_mu, p_sig2, p_rho)
	colnames(RES) <- c("mu", "sig2", "rho")
	RES
}


###data simulation

##correlation matrix
S <- matrix(rho*sig2, q, q)
diag(S) <- sig2
##cholesky
cholS <- chol(S)

##generation from multivariate normal
set.seed(3)
Y <- matrix(rnorm(n*q), n, q)%*%cholS

##compute pairwise score conributions
gr <- pscore_theta(theta, Y)

OBJ.Ib <- Iboot(gradient=gr, B=500, M=500, kB=0.01, alpha=c(0.1, 0.05, 0.01), seed=3)

##critical values for testing H0: theta=(0, 1, 0.5), H1: theta!=(0, 1, 0.5)
OBJ.Ib
summary(OBJ.Ib)

##sampe size too small: convex hull failure at the beginning
n <- 10
set.seed(3)
Y <- matrix(rnorm(n*q), n, q)%*%cholS

##compute pairwise score conributions
gr <- pscore_theta(theta, Y)

OBJ.Ib <- Iboot(gradient=gr, B=500, M=500, kB=0.01, alpha=c(0.1, 0.05, 0.01), seed=3)
OBJ.Ib
summary(OBJ.Ib)
}
}
\keyword{bootstrap}
\keyword{prepivoting}
\keyword{weighted bootstrap}
\keyword{double bootstrap}
\keyword{bootstrap iteration}


