\name{TailPP}
\alias{TailPP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Tail posterior probability for BGmix output.}
\description{
For differential expression models with unstructured priors (no mixture prior), calculates tail posterior probabality and FDR, and plots a histogram. Uses whole posterior distributions of likelihood parameters (found by 'ccTrace') and posterior means of hyperparameters (found by 'ccParams').
}
\usage{
TailPP(res, nreps, params, paired=F, alpha=0.05, N = 5000, prec=0.05, p.cut = 0.7, plots  = T, plot.pi0=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{ list object output from 'ccTrace'}
  \item{nreps}{ vector length 2 containing the number of replicates in each condition}
  \item{params}{ list object output from 'ccParams'}
  \item{paired}{logical. TRUE for paired design, FALSE for unpaired.}
  \item{alpha}{ parameter of the tail posterior probability (1-alpha/2 quantile)}
  \item{N}{ simulation size for tail posterior probability under H0}
  \item{prec}{ parameter used when estimating CDF of tail posterior probability under H0 }
  \item{p.cut}{ calculate FDR only for cutoffs on tail posterior probability > p.cut }
  \item{plots}{ logical. if TRUE, makes plots of the histogram of tail posterior probability with the null density  and of FDR
}
  \item{plot.pi0}{ logical. if TRUE, diagnostic plot of the estimated pi0 at different locations and the median estimate }
}
\value{
  \item{tpp }{vector of tail posterior probabilities with parameter alpha, one per gene}
  \item{FDR }{(smoothed) estimate of FDR for all (distinct) cutoffs > p.cut}
  \item{pi0 }{estimated proportion of observations under the null}
  \item{pp0 }{simulations under the null}
}
\author{ Natalia Bochkina}
\references{Bochkina N., Richardson S. (2007) 
            Tail posterior probability for inference in pairwise and multiclass gene expression data. Biometrics. 
            http://www.blackwell-synergy.com/doi/abs/10.1111/j.1541-0420.2006.00807.x}

\seealso{ \code{\link{FDRplotTailPP}},\code{\link{histTailPP}},\code{\link{EstimatePi0}}}

\examples{

 data(ybar, ss)
 nreps <- c(8,8)

## Note this is a very short MCMC run!
## For good analysis need proper burn-in period.
 outdir <- BGmix(ybar, ss, nreps, jstar=-1, nburn=0, niter=100, nthin=1)

 params <- ccParams(outdir)  
 res <-  ccTrace(outdir)
  
 tpp.res <- TailPP(res, nreps, params, plots  = FALSE)
 histTailPP(tpp.res)
 FDRplotTailPP(tpp.res, plot.TP = TRUE)

}

\keyword{ htest }
