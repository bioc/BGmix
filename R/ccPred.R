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

"ccPred" <-
function(filedir, q.partial=T, q.trace=F, quiet=T){

summ <- ccSummary(filedir)
 
ngenes <- summ$ngenes
nconds <- summ$nconds
ncomps <- summ$ncomps
niter <- summ$niter
nthin <- summ$nthin
  
nsamples <- niter/nthin

pval.ss.post <- matrix(scan(file.path(filedir,"pval_post_ss.txt"),quiet=quiet),ncol=nconds,byrow=T)
pval.ss.mix <- matrix(scan(file.path(filedir,"pval_mix_ss.txt"),quiet=quiet),ncol=nconds,byrow=T)
pval.ybar.post <- matrix(scan(file.path(filedir,"pval_post_ybar.txt"),quiet=quiet),ncol=ncomps,byrow=T)
pval.ybar.mix1 <- matrix(scan(file.path(filedir,"pval_mix1_ybar.txt"),quiet=quiet),ncol=ncomps,byrow=T)
pval.ybar.mix2 <- matrix(scan(file.path(filedir,"pval_mix2_ybar.txt"),quiet=quiet),ncol=ncomps,byrow=T)
pval.ybar.mix3 <- matrix(scan(file.path(filedir,"pval_mix3_ybar.txt"),quiet=quiet),ncol=ncomps,byrow=T)
if(q.partial==T){
  pval.ss.part <- matrix(scan(file.path(filedir,"pval_partial_ss.txt"),quiet=quiet),ncol=nconds,byrow=T)
  pval.ybar.part <- scan(file.path(filedir,"pval_partial_ybar.txt"))
}
else{
  pval.ss.part  <- NULL
  pval.ybar.part  <- NULL
}

if(q.trace){
  ybar.pred1 <- array(scan(file.path(filedir,"trace_ybar_pred1.txt")),dim=c(nconds,ngenes,nsamples))
  print("got ybar.pred1")
  ybar.pred2 <- array(scan(file.path(filedir,"trace_ybar_pred2.txt")),dim=c(nconds,ngenes,nsamples))
  print("got ybar.pred2")
  ybar.pred3 <- array(scan(file.path(filedir,"trace_ybar_pred3.txt")),dim=c(nconds,ngenes,nsamples))
  print("got ybar.pred3")
  ybar.pred4 <- array(scan(file.path(filedir,"trace_ybar_pred4.txt")),dim=c(nconds,ngenes,nsamples))
  print("got ybar.pred4")
  ss.pred1 <- array(scan(file.path(filedir,"trace_ss_pred1.txt")),dim=c(nconds,ngenes,nsamples))
  print("got ss.pred1")
  ss.pred2 <- array(scan(file.path(filedir,"trace_ss_pred2.txt")),dim=c(nconds,ngenes,nsamples))
  print("got ss.pred2")
}
else{
  ybar.pred1 <- NULL
  ybar.pred2 <- NULL
  ybar.pred3 <- NULL
  ybar.pred4 <- NULL
  ss.pred1 <- NULL
  ss.pred2 <- NULL
}

return(list(pval.ss.post=pval.ss.post,pval.ss.mix=pval.ss.mix,
            pval.ybar.post=pval.ybar.post,pval.ybar.mix1=pval.ybar.mix1,
            pval.ybar.mix2=pval.ybar.mix2,pval.ybar.mix3=pval.ybar.mix3,
            pval.ss.part=pval.ss.part,pval.ybar.part=pval.ybar.part, 
            ybar.pred1=ybar.pred1,ybar.pred2=ybar.pred2,ybar.pred3=ybar.pred3,
            ybar.pred4=ybar.pred4,ss.pred1=ss.pred1,ss.pred2=ss.pred2))

}

