\name{histTailPP}
\alias{histTailPP}
\title{Histogram plot for tail posterior probability}
\description{
 Plots a histogram of tail posterior probability with its density under the null hypothesis
}
\usage{
histTailPP(tpp.res, bw=0.05, xlim=c(0,1),nc=10)
}
\arguments{
  \item{tpp.res}{output of TailPP}
  \item{bw}{bandwidth for kernel estimate of the null density}
  \item{xlim}{limits on the x axis}
  \item{nc}{number of bins of the histogram}
}
%\details{}
%\value{}
\references{Bochkina N., Richardson S. (2007) 
            Tail posterior probability for inference in pairwise and multiclass gene expression data. Biometrics. }
\author{Natalia Bochkina}
%\note{}
\seealso{ \code{\link{TailPP}}, \code{\link{FDRplotTailPP}},\code{\link{EstimatePi0}}}


\examples{


data(ybar, ss)
 nreps <- c(8,8)

## Note this is a very short MCMC run!
## For good analysis need proper burn-in period.
 outdir <- BGmix(ybar, ss, nreps, jstar=-1, nburn=0, niter=100, nthin=1)

 params <- ccParams(outdir)  
 res <-  ccTrace(outdir)
  
 tpp.res <- TailPP(res, nreps, params, plots  = FALSE)
 histTailPP(tpp.res, bw=0.04, xlim=c(0,1), nc=10)


}
\keyword{ hplot }
