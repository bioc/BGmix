\name{FDRplotTailPP}
\alias{FDRplotTailPP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot of FDR for tail posterior probability}
\description{
 Plots smoothed FDR vs tail posterior probability or vs the number of differentially expressed (DE) genes
}
\usage{
FDRplotTailPP(tpp.res, nmax = sum(! is.na(tpp.res$FDR)), plot.TP = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tpp.res}{output of TailPP}
  \item{nmax}{ maximum size of the list of DE genes}
  \item{plot.TP}{logical. If TRUE FDR is plotted, otherwise the number of false positives is plotted vs the number of differentially expressed genes}
}
%\details{}
%\value{}
\references{Bochkina N., Richardson S. (2007) 
            Tail posterior probability for inference in pairwise and multiclass gene expression data. Biometrics. }
\author{Natalia Bochkina}
%\note{}
\seealso{ \code{\link{TailPP}}, \code{\link{histTailPP}},\code{\link{EstimatePi0}}}


\examples{


 data(ybar, ss)
 nreps <- c(8,8)

## Note this is a very short MCMC run!
## For good analysis need proper burn-in period.
 outdir <- BGmix(ybar, ss, nreps, jstar=-1, nburn=0, niter=100, nthin=1)

 params <- ccParams(outdir)  
 res <-  ccTrace(outdir)
  
 tpp.res <- TailPP(res, nreps, params, plots  = FALSE)
 FDRplotTailPP(tpp.res, plot.TP = TRUE)


}
\keyword{ hplot }
