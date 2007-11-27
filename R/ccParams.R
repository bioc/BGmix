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

"ccParams" <-
function(filedir, q.beta=T, q.sig=T, q.z=T, quiet=T){

summ <- ccSummary(filedir)
 
ngenes <- summ$ngenes
nconds <- summ$nconds
neffects <- summ$neffects
ncomps <- summ$ncomps
ntau <- summ$ntau
move.choice.aa <- summ$move.choice.aa
move.choice.lam <- summ$move.choice.lam
move.choice.eta <- summ$move.choice.eta
  
### options
q.aa <- move.choice.aa==1
if(is.null(move.choice.eta)) q.eta <- FALSE
else q.eta <- move.choice.eta==1
if(is.null(move.choice.lam)) q.lam <- FALSE
else q.lam <- move.choice.lam==1


if(q.beta){
  mbeta <- matrix(scan(file.path(filedir,"mean_beta.txt"),quiet=quiet),ncol=neffects,byrow=T)
  print("got beta")
}
else mbeta <- NULL

if(q.sig){
  if(ntau==1){
    msig2 <- scan(file.path(filedir,"mean_sig2.txt"),quiet=quiet)
    mtau <- scan(file.path(filedir,"mean_tau.txt"),quiet=quiet)
    print("got sig2")
  }
  else{
    msig2 <- matrix(scan(file.path(filedir,"mean_sig2.txt"),quiet=quiet),ncol=ntau,byrow=T)
    mtau <- matrix(scan(file.path(filedir,"mean_tau.txt"),quiet=quiet),ncol=ntau,byrow=T)
    print("got sig2")
  }
}
else{
  msig2 <- NULL
  mtau <- NULL
}
  
if(q.z){
  mzg <- scan(file.path(filedir,"mean_zg.txt"),quiet=quiet)
  print("got zg")
}
else mzg <- NULL

###################### scalars always read in (if they exist)

mbb <- scan(file.path(filedir,"mean_bb.txt"),quiet=quiet)
if(q.aa) maa <- scan(file.path(filedir,"mean_aa.txt"),quiet=quiet)
else maa <- NULL

if(q.eta) meta <- scan(file.path(filedir,"mean_eta.txt"),quiet=quiet)
else meta <- NULL

if(q.lam) mlambda <- scan(file.path(filedir,"mean_lambda.txt"),quiet=quiet)
else mlambda <- NULL

mwtc <- scan(file.path(filedir,"mean_wtc.txt"),quiet=quiet)

pc <- matrix(scan(file.path(filedir,"prob_class.txt"),quiet=quiet),ncol=ncomps,byrow=T)


return(list(mbeta=mbeta,msig2=msig2,mbb=mbb,maa=maa,mtau=mtau,mwtc=mwtc,
            mzg=mzg,meta=meta,mlambda=mlambda,pc=pc))

}

