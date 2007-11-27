\name{ccPred}
\alias{ccPred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read predictive quantities output from BGmix.}
\description{
Reads predictive p-values from files output from BGmix. Also (optionally) reads posterior predictive distributions of data.
}
\usage{
ccPred(filedir, q.partial = T, q.trace = F, quiet = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filedir}{character. The name of the output directory created by BGmix.}
  \item{q.partial}{ logical. Read partial predictive p-values?}
  \item{q.trace}{ logical. Read posterior predictive distributions of data?}
  \item{quiet}{ logical. Parameter passed to 'scan'. (If false,
'scan' prints details of number of items read in.)}
}
\value{
  \item{pval.ss.post}{matrices no. genes x no. conditions. Posterior
    predictive p-values for sum of squares for each gene in each condition.}
  \item{pval.ss.mix}{matrices no. genes x no. conditions. Mixed predictive p-values for sum of squares for each gene in each condition.}
  \item{pval.ss.part}{matrices no. genes x no. conditions. Partial predictive p-values for sum of squares for each gene in each condition.}
  \item{pval.ybar.post}{matrices no. genes x no. mixture components. Posterior predictive p-values for ybar for each gene in each mixture component.}
  \item{pval.ybar.mix2}{matrices no. genes x no. mixture components. Mixed predictive p-values for ybar for each gene in each mixture component.}
  \item{pval.ybar.part}{matrices no. genes x no. mixture components. Partial predictive p-values for ybar for each gene in each mixture component.}
  \item{ybar.pred1}{Posterior predictive distribution of ybar.}
  \item{ybar.pred3}{Mixed predictive distribution of ybar.}
  \item{ss.pred1}{Posterior predictive distribution of sums of squares.}
  \item{ss.pred2}{Mixed predictive distribution of sums of squares.}
}
\author{ Alex Lewin}
\note{Additional output: pval.ybar.mix1 and pval.ybar.mix3 are alternative versions of mixed predictive p-values (currently not used). Also, ybar.pred2 and ybar.pred4 are the corresponding alternative mixed predictive distributions for ybar.
}
\keyword{ manip }
\examples{
## Note this is a very short MCMC run!
## For good analysis need proper burn-in period.
data(ybar,ss)
outdir <- BGmix(ybar, ss, c(8,8), nburn=0, niter=100, nthin=1)
pred <- ccPred(outdir)
}
