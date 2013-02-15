\name{Iboot-package}
\alias{Iboot-package}
\alias{Ibootpack}
\docType{package}
\title{
Iboot: iterated bootstrap tests and confidence sets}
\description{
	The package provides a computationally efficient and general algorithm to obtain iterated bootstrap tests and confidence sets based on the unstudentised version of the Rao statistic for a \eqn{p}-dimensional parameter. The outer and inner level of resampling required to obtain respectively the simple and the re-calibrated bootstrap critical values (at the null hypothesys) are performed in a weighted fashion. The particular choice of the resampling weights allows to obtain accurate re-calibrated critical values with one level of bootstrap iteration only  (Lee and Young, 2003).

The algorithm is particularly efficient as it combines a deterministic stopping rule (Nankervis, 2005) and a computationally convenient statistic to bootstrap on (Lunardon, 2013).
}
\details{
Function \code{Iboot} is merely an R wrapper to call a set of foreign functions all written in C language so that computational efficiency is increased. Some C routines are borrowed from R sources: numerical optimisation and sorting relies on \code{lbfgsb} and \code{revsort} located in "/src/main/optim.c" and "/src/main/sort.c", respectively. The function \code{ProbSampleReplace} for sampling with unequal probabilities has been slightly modified to cut down the number of unnecessary operations for bootstrap resamplings.
}

\author{
Nicola Lunardon <nicola.lunardon@econ.units.it>.

Maintainer: 
Nicola Lunardon <nicola.lunardon@econ.units.it>.
}
\references{
Lee, S., Young, A. (2003). Prepivoting by weighted bootstrap iteration. \emph{Biometrika}, \bold{90}, 393--410.

Lunardon, N. (2013). Prepivoting composite score statistics by weighted bootstrap iteration. E-print: arXiv/1301.7026.

Nankervis, J. (2005). Computational algorithms for double bootstrap confidence intervals. \emph{Computational statistics & data analysis}, \bold{45}, 461--475.
}

\keyword{ package }
\keyword{ bootstrap }
\keyword{ iterated bootstrap }
\keyword{ double bootstrap }
\keyword{ bootstrap calibration }

\seealso{ 
\code{\link[simpleboot:one.boot]{one.boot}}, \code{\link[boot:boot]{boot}}, \code{\link[stats:stats]{stats}}.
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

####Example 2: variance of a normal with known location
n <- 20
mu <- 1
sig2 <- 1

set.seed(1)
##contributions obtained from the score function
gr <- ( rnorm(n, mu, sqrt(sig2)) - mu )^2/sig2 - 1

OBJ.Ib <- Iboot(gradient=gr, B=500, M=500, kB=0.01, alpha=c(0.1, 0.05, 0.01), seed=3)

##critical values for testing H0: sig2=1, H1: sig2!=1
OBJ.Ib
summary(OBJ.Ib)
}
