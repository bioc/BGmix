\name{BGmix-package}
\alias{BGmix-package}
\docType{package}
\title{
BGmix fits a variety of Bayesian hierarchical models for finding differential gene expression between 2 or more experimental conditions.
}
\description{
BGmix uses a C++ routine to fit the chosen model via an MCMC algorithm. Files are written to a sub-directory in the working directory. The package includes R functions for reading the results into R, and several plotting functions and functions for estimating error rates.
}
\details{
\tabular{ll}{
Package: \tab BGmix\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-02-01\cr
License: \tab GPL\cr
}
See Vignette for details of how to use this package (use openVignette()).
}
\author{
Alex Lewin and Natalia Bochkina

Maintainer: Alex Lewin <a.m.lewin@imperial.co.uk>
}
\references{
  Lewin, A., Bochkina, N. and Richardson, S. (2007), Fully Bayesian
  mixture model for differential gene expression: simulations and model
  checks.
  \url{http://www.bgx.org.uk/publications.html}}
\keyword{ models }
\examples{
## Note this is a very short MCMC run!
## For good analysis need proper burn-in period.
data(ybar,ss)
outdir <- BGmix(ybar, ss, c(8,8), nburn=0, niter=100, nthin=1,trace.pred=1)

## Basic plot of parameters
params <- ccParams(outdir)
plotBasic(params,ybar,ss)

## plots of FDR and related quantities
fdr <- calcFDR(params)
par(mfrow=c(1,2))
plotFDR(fdr)

## plots of Bayesian p-values
## for predictive checks of mixture prior 
pred <- ccPred(outdir,q.trace=TRUE)
plotPredChecks(pred$pval.ybar.mix2,params$pc,probz=0.5)

## plots of predictive density superimposed on data
plotMixDensity(params,pred,ybar,ss)
}
