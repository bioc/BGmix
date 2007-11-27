\name{plotTrace}
\alias{plotTrace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Trace plots for BGmix output.}
\description{
Trace plots are plotted for all scalar parameters. Optionally, traces are plotted for parameters indexed by genes, but for selected genes only.
}
\usage{
plotTrace(res, q.beta = T, q.sig = T, q.z = T, ind.genes = (1:3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{ list object output from 'ccTrace' }
  \item{q.beta}{ logical. Plot trace of beta (gene effect) parameters?}
  \item{q.sig}{ logical. Plot trace of gene variances?}
  \item{q.z}{ logical. Plot trace of gene allocation parameters?}
  \item{ind.genes}{ indices of genes for which to plot gene parameters.}
}
\author{ Alex Lewin}
\keyword{ hplot }
\examples{
## Note this is a very short MCMC run!
## For good analysis need proper burn-in period.
data(ybar,ss)
outdir <- BGmix(ybar, ss, c(8,8), nburn=0, niter=100, nthin=1)
tr <- ccTrace(outdir)
plotTrace(tr)
plotTrace(tr,q.beta=TRUE,q.sig=FALSE,q.z=FALSE,ind.genes=1)
plotTrace(tr,q.beta=FALSE,q.sig=FALSE,q.z=TRUE,ind.genes=sample(1:1000,5))
}
