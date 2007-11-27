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

"ccTrace" <-
function(filedir, q.beta=T, q.sig=T, q.z=T, quiet=T){

summ <- ccSummary(filedir)
  
ngenes <- summ$ngenes
nconds <- summ$nconds
neffects <- summ$neffects
ncomps <- summ$ncomps
ntau <- summ$ntau
niter <- summ$niter
nthin <- summ$nthin
move.choice.aa <- summ$move.choice.aa
move.choice.lam <- summ$move.choice.lam
move.choice.eta <- summ$move.choice.eta
  
nsamples <- niter/nthin

### options
q.aa <- move.choice.aa==1
if(is.null(move.choice.eta)) q.eta <- FALSE
else q.eta <- move.choice.eta==1
if(is.null(move.choice.lam)) q.lam <- FALSE
else q.lam <- move.choice.lam==1

##################################################
## beta is factors x genes x samples
## sig2 is ntau x genes x samples
## bb is samples x conds
##################################################

################### options for reading large matrices
  
if(q.beta){
  beta <- array(scan(file.path(filedir,"trace_beta.txt"),quiet=quiet),dim=c(neffects,ngenes,nsamples))
  print("got beta")
}
else beta <- NULL

if(q.sig){
  if(ntau==1){
    sig2 <- array(scan(file.path(filedir,"trace_sig2.txt"),quiet=quiet),dim=c(ngenes,nsamples))
    print("got sig2")
  }
  else{
    sig2 <- array(scan(file.path(filedir,"trace_sig2.txt"),quiet=quiet),dim=c(ntau,ngenes,nsamples))
    print("got sig2")
  }
}
else sig2 <- NULL

if(q.z){
  zg <- matrix(scan(file.path(filedir,"trace_zg.txt"),quiet=quiet),ncol=ngenes,byrow=T)
  print("got zg")
}
else zg <- NULL

###################### scalars always read in (if they exist)

if(ntau==1){
  bb <- scan(file.path(filedir,"trace_bb.txt"),quiet=quiet)
  print("got bb")
  if(q.aa){ 
    aa <- scan(file.path(filedir,"trace_aa.txt"),quiet=quiet)
    print("got aa")
  }
  else aa <- NULL
}
else{
  bb <- matrix(scan(file.path(filedir,"trace_bb.txt"),quiet=quiet),ncol=ntau,byrow=T)
  print("got bb")
  if(q.aa){ 
    aa <- matrix(scan(file.path(filedir,"trace_aa.txt"),quiet=quiet),ncol=ntau,byrow=T)
    print("got aa")
  }
  else aa <- NULL
}

if(q.eta){ 
  eta <- array(scan(file.path(filedir,"trace_eta.txt"),quiet=quiet),dim=c(ncomps-1,nsamples))
  print("got eta")
}
else eta <- NULL

if(q.lam){ 
  lambda <- array(scan(file.path(filedir,"trace_lambda.txt"),quiet=quiet),dim=c(ncomps-1,nsamples))
  print("got lambda")
}
else lambda <- NULL

wtc <- matrix(scan(file.path(filedir,"trace_wtc.txt"),quiet=quiet),ncol=ncomps,byrow=T)
print("got wtc")


return(list(summ=summ,eta=eta,lambda=lambda,aa=aa,bb=bb,wtc=wtc,beta=beta,sig2=sig2,zg=zg))

}

