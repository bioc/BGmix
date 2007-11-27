#  This file is part of BGmix, a fully Bayesian model for
#  differential expression.
#  Copyright 2007 Alex Lewin <a.m.lewin@imperial.ac.uk>
#
#  BGmix is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License, version 2, as
#  published by the Free Software Foundation.
#
#  BGmix is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

"BGmix" <-
function(ybar, ss, nreps, neffects=2, xx=matrix(c(1,1,-0.5,0.5),ncol=2,byrow=T),
         ntau=NULL, indtau=NULL, jstar=1, niter=10000, nburn=10000, nthin=10,  
	seed=12345, move.choice.bz=4, move.choice.aa=1, move.choice.lam=0,
	    move.choice.tau=1, move.choice.eta=1, trace.out=1, trace.pred=0,
	 sig.aa=0.1, tau.eps=50,datafilename.ybar=NULL, xfilename=NULL,
         itfilename=NULL, rundir="."){

############ variables which are set automatically from others

if(sum(nreps==0)>0) stop("Sorry, this program only works with data coming from >1 replicates")
  
ngenes <- dim(ybar)[1]
nconds <- dim(ybar)[2]
nrepstot <- sum(nreps)
ydata <- rep(0,ngenes*nrepstot)

if(is.null(datafilename.ybar)) datafilename.ybar <- substitute(ybar)
if(is.null(xfilename)) xfilename <- "x def. in R"
if(is.null(itfilename)) itfilename <- "indtau def. in R"

if(is.null(ntau)) ntau <- nconds
if(is.null(indtau)) indtau <- 0:(nconds-1)

############### variables which I've taken out to make the list shorter

aa.init=0.8; gg=0.01; hh=0.01; tau.var=0.001; eta.up.init=3; eta.down.init=3 
lambda.up.init=1.5; lambda.down.init=1.5; aa.eta=1; bb.eta=1
lam1=0.2; lam2=5; nlam=25; nu1=1; nu2=1; nu0=1
zg.init=0; beta.init=0.1; tau.init=1.1; bb.init=0.2

############### variables which I don't want people to change
  
move.choice.cut <- 1
like.choice <- 1
df <- 4
ncomps <- 3

######################## print out and check options

if(dim(xx)[1]!=neffects) stop("neffects, xx inconsistent, xx should be neffects x nconds")
if(dim(xx)[2]!=nconds) stop("nconds, xx inconsistent, xx should be neffects x nconds")

   
if(nburn>0 & nburn<100) stop("nburn must be >=100 (or 0)")
if(niter>0 & niter<100) stop("niter must be >=100 (or 0)")

if(jstar==-1) print("Flat model (no mixture prior)")
else if(jstar>=0 & jstar<=(ncomps-1)) print(paste("Mixture prior on comp.",jstar+1))
else stop("invalid jstar")

if(jstar!=-1){
  if(move.choice.bz==1){
    print("delta ~ Unif, Gibbs")
    if(move.choice.eta==1) print("eta (limit of Unif) updated. Not very identifiable!")
  }
  else if(move.choice.bz==2){
    print("delta ~ Unif, MH")
    if(move.choice.eta==1) print("eta (limit of Unif) updated. Not very identifiable!")
  }
  else if(move.choice.bz==3 | move.choice.bz==4){
    print("delta ~ Gamma, MH")
    if(move.choice.eta==1) print("eta (scale of Gamma) updated")
    else print("eta (scale of Gamma) NOT updated")
    if(move.choice.lam==1) print("lambda (shape of Gamma) updated. Not always identifiable!")
    else print("lambda (shape of Gamma) not updated")
  }  
  else if(move.choice.bz==5){
    print("delta ~ Gamma + nugget null, MH")
    if(move.choice.eta==1) print("eta (scale of Gamma) updated")
    else print("eta (scale of Gamma) NOT updated")
    if(move.choice.lam==1) print("lambda (shape of Gamma) updated. Not always identifiable!")
    else print("lambda (shape of Gamma) not updated")
  }
  else stop("Invalid move.choice.bz (should be 1,...,5)")
}
  
if(like.choice==1) print("Normal Likelihood")
else if(like.choice==2) print("T Likelihood")
else stop("Invalid like.choice (should be 1 or 2)")

if(move.choice.tau==1){
  print("tau ~ Gamma")
  if(move.choice.cut!=1) print("WARNING: cut fn used for tau: This facility not maintained.")
  if(move.choice.aa==1) print("a (prior for tau) is updated")
  else print("a (prior for tau) NOT updated")
}
else{
  print("tau ~ logNormal")
  print("a (prior for tau) is updated")
}

if(trace.out==1) print("trace output for parameters")
else print("no trace for parameters")

if(trace.pred==1) print("trace output for predicted data")
else print("no trace for predicted data")

dirname <- fitBGmix(ngenes, nconds, neffects, ncomps, ntau, niter, nburn, nthin, jstar, 
	seed, zg.init, nlam, move.choice.bz, move.choice.cut, move.choice.aa, move.choice.lam,
	    move.choice.tau, move.choice.eta, like.choice, trace.out, trace.pred,
	aa.init, gg, hh, tau.var, eta.up.init, eta.down.init, 
	lambda.up.init, lambda.down.init, aa.eta, bb.eta, lam1, lam2, nu1, nu2, 
	nu0, df, beta.init, tau.init, bb.init, sig.aa, tau.eps, ydata,
	ybar, ss, nreps, nrepstot, xx, indtau, datafilename.ybar, xfilename, itfilename, rundir)
  
print(paste("Output directory is",dirname))

return(dirname)

}

