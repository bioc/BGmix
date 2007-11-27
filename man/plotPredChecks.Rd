\name{plotPredChecks}
\alias{plotPredChecks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plots of predictive checks for mixture prior.}
\description{
Histograms and q-q plots of predictive p-values for the mixture prior. Separate plots are given for each mixture component, using only genes with high posterior probability of being classified into the relevant component. 
}
\usage{
plotPredChecks(pvals, pc, probz = 0.8, label = "", breaks = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvals}{ matrix of predictive p-values output by 'ccPred' (NB, not the whole list object, just the matrix of p-values)}
  \item{pc}{ matrix of posterior classification probabilities output by 'ccParams' (NB, not the whole list object, just the matrix of probabilities)}
  \item{probz}{ threshold on posterior probabilities for including genes in each mixture component plot}
  \item{label}{ title used on histograms}
  \item{breaks}{ argument input to histogram}
}
\author{ Alex Lewin}
\keyword{ hplot }
\examples{
## Note this is a very short MCMC run!
## For good analysis need proper burn-in period.
data(ybar,ss)
outdir <- BGmix(ybar, ss, c(8,8), nburn=0, niter=100, nthin=1)
params <- ccParams(outdir)
pred <- ccPred(outdir)
plotPredChecks(pred$pval.ybar.mix2,params$pc,probz=0.5)
}
