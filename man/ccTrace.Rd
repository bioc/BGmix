\name{ccTrace}
\alias{ccTrace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read trace files from BGmix}
\description{
Reads output files containing whole posterior distributions from BGmix. Also calls 'ccSummary', and outputs model options.
}
\usage{
ccTrace(filedir, q.beta = T, q.sig = T, q.z = T, quiet = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filedir}{character. The name of the output directory created by BGmix.}
  \item{q.beta}{ logical. Read beta values?}
  \item{q.sig}{ logical. Read gene variances?}
  \item{q.z}{ logical. Read z values?}
  \item{quiet}{ logical. Parameter passed to 'scan'. (If false,
    'scan' prints details of number of items read in.)}
  }
\value{
  \item{summ }{list object output by 'ccSummary'}
  \item{eta }{matrix (no. components -1) x no. MCMC samples. Posterior of mixture component parameters (eta).}
  \item{lambda }{matrix (no. components -1) x no. MCMC samples. Posterior of mixture component parameters (lambda).}
  \item{aa }{matrix no. MCMC samples x no. variances. Posterior of variance hyperparameters (a).}
  \item{bb }{matrix no. MCMC samples x no. variances. Posterior of variance hyperparameters (b).}
  \item{wtc }{matrix no. MCMC samples x no. mixture components. Posterior of mixture weights.}
  \item{beta }{matrix no. effects x no. genes x no. MCMC samples. Posterior of gene effects.}
  \item{sig2 }{matrix no. variances x no. genes x no. MCMC samples. Posterior of gene variances.}
  \item{zg }{matrix no. MCMC samples x no. genes. Posterior of gene allocations.}
}
\author{ Alex Lewin}
\keyword{ manip }
\examples{
## Note this is a very short MCMC run!
## For good analysis need proper burn-in period.
data(ybar,ss)
outdir <- BGmix(ybar, ss, c(8,8), nburn=0, niter=100, nthin=1)
tr <- ccTrace(outdir)
}
