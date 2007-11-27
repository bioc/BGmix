\name{plotBasic}
\alias{plotBasic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Basic plots of BGmix parameters and data.}
\description{
Plots gene effects and variances versus their corresponding data sufficient statistics (to show the effect of smoothing and shrinkage). Also plots "volcano plots": posterior probabilities of being classified in  each mixture component versus the log fold change parameters.
}
\usage{
plotBasic(res, ybar, ss, q.mean = T, q.diff = T, q.sig = T, q.volcano = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{ list object output from 'ccParams'}
  \item{ybar}{ ybar data (see \code{\link{BGmix}} help for details)}
  \item{ss}{ ss data (see \code{\link{BGmix}} help for details)}
  \item{q.mean}{ logical. Include mean plot?}
  \item{q.diff}{ logical. Include log fold change plot?}
  \item{q.sig}{ logical. Include variance plot?}
  \item{q.volcano}{ logical. Include volcana plot (posterior classification v. fold change)? }
}
\details{
Note this plotting function is designed for model output from the unpaired differential expression design.}
\value{
No value is returned to \R. Results from BGmix model are output to files. 
}
\author{ Alex Lewin}
\keyword{ hplot }
\examples{
## Note this is a very short MCMC run!
## For good analysis need proper burn-in period.
data(ybar,ss)
outdir <- BGmix(ybar, ss, c(8,8), nburn=0, niter=100, nthin=1)
params <- ccParams(outdir)
plotBasic(params,ybar,ss)
}
