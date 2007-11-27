\name{ccParams}
\alias{ccParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read posterior means and classification probabilities from BGmix}
\description{
Reads output files containing posterior means from BGmix AND reads posterior probabilities of each gene being classified in the null mixture component.
}
\usage{
ccParams(filedir, q.beta = T, q.sig = T, q.z = T, quiet = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filedir}{character. The name of the output directory created by BGmix.}
  \item{q.beta}{ logical. Read beta values?}
  \item{q.sig}{ logical. Read gene variance parameters?}
  \item{q.z}{ logical. Read z values?}
  \item{quiet}{ logical. Parameter passed to 'scan'. (If false, 'scan' prints details of number of items read in.)}
}
\value{
  \item{mbeta }{matrix no. genes x no. effects. Posterior means of gene effect parameters (usually gene means and log fold changes).}
  \item{msig2}{matrix no. genes x no. variances. Posterior means of gene variances.}
  \item{mbb}{vector of hyperparameters (b) for gene variances (posterior means).}
  \item{maa}{vector of hyperparameters (a) for gene variances (posterior means).}
  \item{mtau}{matrix no. genes x no. conditions. Posterior means of gene precisions.}
  \item{mwtc}{vector of posterior mean mixture weights}
  \item{mzg}{vector of posterior mean allocation for each gene}
  \item{meta}{vector of mixture parameters (eta)}
  \item{mlambda}{vector of mixture parameters (lambda)}
  \item{pc}{matrix no. genes x no. mixture components. Posterior probability for each gene of being classified into each mixture component.}
}
\author{Alex Lewin}
\keyword{ manip }
\examples{
## Note this is a very short MCMC run!
## For good analysis need proper burn-in period.
data(ybar,ss)
outdir <- BGmix(ybar, ss, c(8,8), nburn=0, niter=100, nthin=1)
params <- ccParams(outdir)
}

