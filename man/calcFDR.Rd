\name{calcFDR}
\alias{calcFDR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimate the FDR (false discovery rate) and related quantities for BGmix output.}
\description{
Given a threshold on the posterior probabilities, genes are declared as null or differentially expressed. For any given threshold, the FDR  (false discovery rate) and FNR (false non-discovery rate) can be estimated using the posterior probabilities. Estimated numbers of false positives and false negatives are also output.
}
\usage{
calcFDR(res, pcut = seq(0.01,0.5,0.01), true.z = NULL, q.print = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{ list object output from ccParams (this includes the posterior classification probabilities) }
  \item{pcut}{ scalar or vector of thresholds for which to estimate FDR etc. }
  \item{true.z}{vector of true classifications (if known, eg. for simulated data) }
  \item{q.print}{ Print FDR etc. when pcut is a vector?}
}
\details{
  If the true classification is known, it can be given as true.z, and
  the true FDR etc. for the threshold probability can be calculated.
}
\value{
  \item{fdr.est, fnr.est}{scalars or vectors of estimated FDR, FNR}
  \item{fp.est, fn.est}{scalars or vectors of estimated no. false positives,
  no. false negatives}
  \item{fdr.true, fnr.true}{scalars or vectors of true FDR, FNR}
  \item{fp.true, fn.true}{scalars or vectors of true no. false positives,
  no. false negatives}
  \item{npos, nneg}{scalars or vectors of no. declared positives,
    no. declared negatives}
  \item{prob.class}{posterior classification probabilites (from the 'res' object input to this function)}
  \item{true.z}{argument to function is output}
  \item{pcut}{argument to function is output}
}
\author{ Alex Lewin}
\keyword{ htest }
\examples{
## Note this is a very short MCMC run!
## For good analysis need proper burn-in period.
data(ybar,ss)
outdir <- BGmix(ybar, ss, c(8,8), nburn=0, niter=100, nthin=1)
params <- ccParams(outdir)
fdr <- calcFDR(params)
}
