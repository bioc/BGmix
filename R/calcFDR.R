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

"calcFDR" <-
function(res,pcut=seq(0.01,0.5,0.01),true.z=NULL,q.print=F){

prob.class <- res$pc[,2]

### prob.class must be a vector of posterior probs of
### being in the null
### true.z must be a vector of 1 for DE, 0 for null

### pcut may be scalar or vector
### output fdr.est etc. are vectors if pcut is vector
  
n <- length(prob.class)
nvec <- length(pcut)

if(nvec==1){
 est.z <- prob.class < pcut
 npos <- sum(est.z==1)
 nneg <- sum(est.z==0)
 if(!is.null(true.z)){
  fp.true <- sum((est.z-true.z)==1)
  fn.true <- sum((true.z-est.z)==1)
  fdr.true <- fp.true/npos
  fnr.true <- fn.true/nneg
  print(paste("True FDR = ",fdr.true))
  print(paste("True FNR = ",fnr.true))
 }
 fp.est <- sum(prob.class*est.z)
 fn.est <- sum((1-prob.class)*(1-est.z))
 fdr.est <- fp.est/npos
 fnr.est <- fn.est/nneg
 print(paste("Estimated FDR = ",fdr.est))
 print(paste("Estimated FNR = ",fnr.est))
}
else{
 ## replicate prob.class nvec times to cf with pcut (length nvec)
 est.p <- t( matrix(rep(prob.class,nvec),ncol=nvec) ) 
 est.z <- est.p < pcut
 npos <- apply( est.z==1 ,1,sum)
 nneg <- apply( est.z==0 ,1,sum) 
 if(!is.null(true.z)){
  fp.true <- apply( (est.z-true.z)==1 ,1,sum)
  fn.true <- apply( (true.z-est.z)==1 ,1,sum)
  fdr.true <- fp.true/npos
  fnr.true <- fn.true/nneg
  if(q.print){
    print(paste("True FDR = ",fdr.true))
    print(paste("True FNR = ",fnr.true))
  }
 }
 fp.est <- apply( est.p*est.z ,1,sum)
 fn.est <- apply( (1-est.p)*(1-est.z) ,1,sum)
 fdr.est <- fp.est/npos
 fnr.est <- fn.est/nneg
 if(q.print){
   print(paste("Estimated FDR = ",fdr.est))
   print(paste("Estimated FNR = ",fnr.est))
 }
}

if(is.null(true.z)){
 fdr.true <- NULL
 fnr.true <- NULL
 fp.true <- NULL
 fn.true <- NULL
}

return(list(fdr.est=fdr.est,fnr.est=fnr.est,fp.est=fp.est,fn.est=fn.est,
            fdr.true=fdr.true,fnr.true=fnr.true,fp.true=fp.true,fn.true=fn.true,
            npos=npos,nneg=nneg,prob.class=prob.class,true.z=true.z,pcut=pcut))

}

