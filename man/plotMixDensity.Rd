\name{plotMixDensity}
\alias{plotMixDensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot predictive density of data.}
\description{
Plot predictive density of data superimposed on histograms of observed
data. Separate plots for ybar and sums of squares.}
\usage{
plotMixDensity(res, predres, ybar, ss)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{ list object output from 'ccParams'}
  \item{predres}{ list object output from 'ccPred' (need q.trace=T in 'ccPred')}
  \item{ybar}{ ybar data (see \code{\link{BGmix}} help for details)}
  \item{ss}{ ss data (see \code{\link{BGmix}} help for details)}
}
\details{
Note that this function is written for the unpaired differential
expression design.
  }
\author{ Alex Lewin}
\keyword{ hplot }
\examples{
## Note this is a very short MCMC run!
## For good analysis need proper burn-in period.
data(ybar,ss)
outdir <- BGmix(ybar,ss,c(8,8),niter=100,nburn=0,nthin=1,trace.pred=1)
pred <- ccPred(outdir,q.trace=TRUE)
params <- ccParams(outdir)
plotMixDensity(params,pred,ybar,ss)
}
