\name{ccSummary}
\alias{ccSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read summary of BGmix output}
\description{
Reads the summary.txt file output by BGmix, containing information about data sets used and model options. This function is called by ccParams, ccTrace and ccPred, therefore users will not in general need to call it directly. 
}
\usage{
ccSummary(filedir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filedir}{character. The name of the output directory created by BGmix.}
}
\value{
A list of scalar values, as follows:
  \item{ngenes,nconds,neffects,ncomps,ntau}{nos. genes, conditions, effects, mixture components, gene variances}
  \item{jstar }{label of effect with mixture prior (labels start at 0)}
  \item{move.choice.bz,move.choice.cut,move.choice.aa,
       move.choice.eta,move.choice.lam,move.choice.tau,
       like.choice,trace.out,trace.pred}{model choice options (see \code{\link{BGmix}} help for details}
  \item{lambda.up.init,lambda.down.init,
            eta.up.init,eta.down.init
}{inital values for eta and lambda (parameters of mixture components)}
}
\author{Alex Lewin}
\keyword{ manip }
\examples{
## Note this is a very short MCMC run!
## For good analysis need proper burn-in period.
data(ybar,ss)
outdir <- BGmix(ybar, ss, c(8,8), nburn=0, niter=100, nthin=1)
summ <- ccSummary(outdir)
}
