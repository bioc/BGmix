\name{FDRforTailPP}
\alias{FDRforTailPP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{FDR for tail posterior probability}
\description{
Calculate the false discovery rate (FDR) for the tail posterior probability
}
\usage{
FDRforTailPP(tpp, a1, a2 = NULL, n.rep1, n.rep2 = NULL, prec = 0.05, p.cut = 0.7, N = 10000, pp0=NULL, plot = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tpp}{vector of tail posterior probabilities }
  \item{a1}{  posterior mean of the shape parameter of the inverse gamma distribution - prior for the variance in condition 1 }
  \item{a2}{  posterior mean of the shape parameter of the inverse gamma distribution - prior for the variance in condition 2 }
  \item{n.rep1}{ number of replicates in condition 1 }
  \item{n.rep2}{ number of replicates in condition 2 }
  \item{prec}{ precision of the estimate of the cumulative distribution function of tail posterior probability under H0 (at points 1 - k*prec, k =1,2,..) }
  \item{p.cut}{ to save time, calculate FDR only for cutoffs on tail posterior probability > p.cut }
  \item{N}{ simulation size for tail posterior probability under H0 }
  \item{pp0}{ a vector of simulated tail posterior probabilities under H0 }
  \item{plot}{ if True, the estimated pi0 at different locations and the median estimate is plotted}
}
%\details{}
\value{
  \item{pi0}{estimate of pi0 - proportion of non-differentially expressed genes}
  \item{ FDR }{estimate of FDR for all (distinct) cutoffs > p.cut}
}
\references{Bochkina N., Richardson S. (2007) 
            Tail posterior probability for inference in pairwise and multiclass gene expression data. Biometrics. 
 }
\author{Natalia Bochkina}
%\note{ }
\seealso{ \code{\link{TailPP}}, \code{\link{FDRplotTailPP}},\code{\link{histTailPP}},\code{\link{EstimatePi0}}}

\examples{


 data(ybar, ss)
 nreps <- c(8,8)

## Note this is a very short MCMC run!
## For good analysis need proper burn-in period.
 outdir <- BGmix(ybar, ss, nreps, jstar=-1, nburn=0, niter=100, nthin=1)

 params <- ccParams(outdir)  
 res <-  ccTrace(outdir)
  
 tpp.res <- TailPP(res, nreps, params, plots  = FALSE)
 FDR.res = FDRforTailPP(tpp.res$tpp, a1 = params$maa[1],
a2 = params$maa[2], n.rep1=nreps[1], n.rep2=nreps[2], p.cut = 0.8)


}
\keyword{ htest }
