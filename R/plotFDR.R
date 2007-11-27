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

"plotFDR" <-
function(res,ylim=NULL,q.plotfnr=F,q.plotpcut=T,q.plotnpos=T,...){

if(is.null(ylim)) ylim <- c(0,0.5)

##### find gene lists for certain FDR values
np1 <- sum(res$fdr.est<=0.05)
np2 <- sum(res$fdr.est<=0.1)

#if(names(dev.cur())=="X11"){
#  X11()
#  par(mfrow=c(2,1))
#}

################################# plot v. pcut

if(q.plotpcut){

###### est. fdr, fnr v. pcut
plot(res$pcut,res$fdr.est,ylim=ylim,type="l",
     xlab="Cut-off on posterior probability",ylab="Estimated FDR",...)
abline(v=res$pcut[np1],h=res$fdr.est[np1],col="red")
abline(v=res$pcut[np2],h=res$fdr.est[np2],col="blue")
xx <- res$pcut[1] + 0.1*(res$pcut[length(res$pcut)]-res$pcut[1])
yy <- ylim[1] + 0.9*(ylim[2]-ylim[1])
legend(xx,yy,leg=c("FDR=5%","FDR=10%"),col=c("red","blue"),lty=c(1,1),bg="white")

if(q.plotfnr) plot(res$pcut,res$fnr.est,ylim=ylim,type="l",
     xlab="Cut-off on posterior probability",ylab="Estimated FNR",...)

###### true fdr, fnr v. pcut
if(!is.null(res$fdr.true)) plot(res$pcut,res$fdr.true,ylim=ylim,type="l",
     xlab="Cut-off on posterior probability",ylab="True FDR",...)
if(!is.null(res$fnr.true)) plot(res$pcut,res$fnr.true,ylim=ylim,type="l",
     xlab="Cut-off on posterior probability",ylab="True FNR",...)

}

################################### plot v. npos

if(q.plotnpos){
  
###### est. fdr, fnr v. no. DE genes
plot(res$npos,res$fdr.est,ylim=ylim,type="l",
     xlab="No. DE genes",ylab="Estimated FDR",...)
abline(v=res$npos[np1],h=res$fdr.est[np1],col="red")
abline(v=res$npos[np2],h=res$fdr.est[np2],col="blue")
xx <- res$npos[1] + 0.1*(res$npos[length(res$npos)]-res$npos[1])
yy <- ylim[1] + 0.9*(ylim[2]-ylim[1])
legend(xx,yy,leg=c("FDR=5%","FDR=10%"),col=c("red","blue"),lty=c(1,1),bg="white")

if(q.plotfnr) plot(res$npos,res$fnr.est,ylim=ylim,type="l",
     xlab="No. DE genes",ylab="Estimated FNR",...)

}

print(paste("No. DE genes for FDR<=5% is ",res$npos[np1]))
print(paste("No. DE genes for FDR<=10% is ",res$npos[np2]))


}

