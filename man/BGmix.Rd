\name{BGmix}
\alias{BGmix}
\title{Fit the BGmix differential expression model. }
\description{
This is the main function of the BGmix package. It calls the C++ code which performs the MCMC to fit the BGmix model.
}
\usage{
BGmix(ybar, ss, nreps, neffects = 2, xx = matrix(c(1, 1, -0.5, 0.5),
ncol = 2, byrow = T), ntau = NULL, indtau = NULL, jstar = 1, niter =
10000, nburn = 10000, nthin = 10, seed = 12345, move.choice.bz = 4,
move.choice.aa = 1, move.choice.lam = 0, move.choice.tau = 1,
move.choice.eta = 1, trace.out = 1, trace.pred = 0, sig.aa = 0.1,
tau.eps = 50, lambda.up.init=1.5, lambda.down.init=1.5,
datafilename.ybar = NULL, xfilename = NULL, itfilename =
NULL, rundir=".")
}
\arguments{
  \item{ybar}{ matrix no. genes x no. experimental conditions. Mean log gene expression for each gene in each condition.}
  \item{ss}{ matrix no. genes x no. experimental conditions. Sample variance of log gene expression for each gene in each condition.}
  \item{nreps}{ vector containing the number of replicate arrays in each experimental condition}
  \item{neffects}{ number of effect parameters per gene (eg. 2 for unpaired differential expression)}
  \item{xx}{ design matrix: no. effects x no. experimental conditions. See Vignette for specification of design matrix. Default is for unpaired differential expression.}
  \item{ntau}{ number of variances per gene }
  \item{indtau}{ label for each condition indicating which variance grouping that condition belongs to. See Vignette for more detail. }
  \item{jstar}{ Label of the effect parameter which has the mixture prior. Labels start at 0, as in C++. If no parameter has a mixture prior, set jstar=-1.}
  \item{niter}{ no. MCMC iterations after burn-in. This must be at least
  100 for the function to work (or else set to zero).}
  \item{nburn}{ no. MCMC iterations for burn-in. This must be at least
    100 for the function to work (or else set to zero).}
  \item{nthin}{ thinning parameter for MCMC iterations}
  \item{seed}{ initial value for random seed}
  \item{move.choice.bz}{ indicates choice of mixture prior: 1 for point mass null + Uniform alternatives, 4 for point mass null + Gamma alternatives, 5 for small Normal null + Gamma alternatives}
  \item{move.choice.aa}{if this is 1, hyperparameter a for gene variances is updated, if this is 0 it is fixed. }
  \item{move.choice.lam}{ if this is 1, hyperparameter lambda for mixture prior is updated, if this is 0 it is fixed. }
  \item{move.choice.tau}{ indicates choice of prior on gene variances: 1 for Inverse Gamma, 2 for log Normal.}
  \item{move.choice.eta}{ if this is 1, hyperparameter eta for mixture prior is updated, if this is 0 it is fixed. }
  \item{trace.out}{if this is 1, output trace of model parameters, if this is 0, no output.}
  \item{trace.pred}{ if this is 1, output trace of predictive quantities, if this is 0, no output.}
  \item{sig.aa}{ step-size in random walk update for a (hyperparameter for gene variances distribution)}
  \item{tau.eps}{Value of epsilon used in the small Normal null mixture
    component.}
  \item{lambda.up.init}{init or fixed value of lambda+ (parameter of
    Gamma mixture component)}
  \item{lambda.down.init}{init or fixed value of lambda- (parameter of
    Gamma mixture component)}
  \item{datafilename.ybar}{ character. Name describing the data set (by default this is taken from the name of the ybar argument).}
  \item{xfilename}{ character. Name describing the design matrix.}
  \item{itfilename}{ character. Name describing the indtau parameter.}
  \item{rundir}{ character. Path for saving output files. A new
    sub-directory is created in the \code{rundir} directory.}
}
\details{
  The C++ code writes a count down on the screen, to give an indication
  of how long the code has to run. Output is written to a sub-directory
  of the working directory. This sub-directory is created automatically,
  and its name is printed by the C++ code to the screen.
}
\value{
  The output directory is returned (character).
}
\references{
  Lewin, A., Bochkina, N. and Richardson, S. (2007), Fully Bayesian
  mixture model for differential gene expression: simulations and model
  checks.
  \url{http://www.bgx.org.uk/publications.html}}
\author{ Alex Lewin}
\keyword{ models }
\examples{
## Note this is a very short MCMC run!
## For good analysis need proper burn-in period.
data(ybar,ss)
BGmix(ybar, ss, c(8,8), nburn=0, niter=1000, nthin=1)
}

