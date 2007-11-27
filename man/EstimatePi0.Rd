\name{EstimatePi0}
\alias{EstimatePi0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Proportion of the variables under the null hypothesis}
\description{
Estimate of the proportion of the variables under the null hypothesis using tail posterior probabilities
}
\usage{
EstimatePi0(tpp, pp0, plot = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tpp}{ observed tail posterior probability }
  \item{pp0}{ a vector of tail posterior probability under H0}
  \item{plot}{ if True, estimated pi0 at different locations and the median estimate is plotted }
}
\details{
  Use Storey (2002) approach to estimate pi0
}
\value{
 estimate of pi0 = proportion of non-differentially expressed genes
}
\references{Bochkina N., Richardson S. (2007) 
            Tail posterior probability for inference in pairwise and multiclass gene expression data. Biometrics (in press). 
         }
\author{Natalia Bochkina}
\seealso{ \code{\link{TailPP}}, \code{\link{FDRplotTailPP}},\code{\link{histTailPP}}}

\examples{

 data(ybar, ss)
 nreps <- c(8,8)

## Note this is a very short MCMC run!
## For good analysis need proper burn-in period.
 outdir <- BGmix(ybar, ss, nreps, jstar=-1, nburn=0, niter=100, nthin=1)

 params <- ccParams(outdir)  
 res <-  ccTrace(outdir)
  
 tpp.res <- TailPP(res, nreps, params, plots  = FALSE)
 pi0 <- EstimatePi0(tpp.res$tpp, tpp.res$pp0)

}

\keyword{ htest }
