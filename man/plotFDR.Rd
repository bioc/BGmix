\name{plotFDR}
\alias{plotFDR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot estimated FDR etc. for BGmix output.}
\description{
Given a threshold on the posterior probabilities, genes are declared as
null or differentially expressed. For any given threshold, the FDR
(false discovery rate) and FNR (false non-discovery rate) can be
estimated using the posterior probabilities. This function plots these
quantities twice, once versus the threshold probabilities, and once versus the number of declared positives. 
}
\usage{
plotFDR(res, ylim = NULL, q.plotfnr = F, q.plotpcut = T, q.plotnpos = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{ list object output from 'calcFDR'}
  \item{ylim}{ optional argument specifying limit for y-axis}
  \item{q.plotfnr}{ Include FNR in plots?}
  \item{q.plotpcut}{ Include the plot of error rates v. threshold on posterior probabilites?}
  \item{q.plotnpos}{ Include the plot of error rates v. no. positives.}
  \item{\dots}{arguments passed to 'plot'}
}
\value{
No value is returned to \R. Results from BGmix model are output to files. }
\author{ Alex Lewin}
\keyword{ hplot }
\examples{
## Note this is a very short MCMC run!
## For good analysis need proper burn-in period.
data(ybar,ss)
outdir <- BGmix(ybar, ss, c(8,8), nburn=0, niter=100, nthin=1)
params <- ccParams(outdir)
fdr <- calcFDR(params)
par(mfrow=c(1,2))
plotFDR(fdr)
}
