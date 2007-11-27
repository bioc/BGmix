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

"plotPredChecks" <-
function(pvals,pc,probz=0.8,label="",breaks=20){
  
xlim <- c(0,1)
  
pval1 <- pvals[,1]
pval2 <- pvals[,2]
pval3 <- pvals[,3]
pc1 <- pc[,1]
pc2 <- pc[,2]
pc3 <- pc[,3]

ind1 <- pc1 >= probz
ind2 <- pc2 >= probz
ind3 <- pc3 >= probz

min.lhs <- min(pval1[ind1])
max.rhs <- max(pval3[ind3])

n1 <- sum(!is.nan(pval1[ind1]))
n2 <- sum(!is.nan(pval2[ind2]))
n3 <- sum(!is.nan(pval3[ind3]))

#print(n1)
#print(n2)
#print(n3)

########### histograms conditional on z, only genes where strong inference on zg

if(sum(ind1)>1){
  h <- hist(pval1[ind1],main=paste(label," z=-1"),xlab="mixed predictive p-value",xlim=xlim,breaks=breaks)
  #mtext(sum(h$counts))
  #title(sub=paste("only genes prob(z)>",probz))
}
else plot(c(0,1),c(0,1),type="n")

if(sum(ind2)>1){
  h <- hist(pval2[ind2],main=paste(label," z=0"),xlab="mixed predictive p-value",xlim=xlim,breaks=breaks)  
  #mtext(sum(h$counts))
  #title(sub=paste("only genes prob(z)>",probz))
}
else plot(c(0,1),c(0,1),type="n")

if(sum(ind3)>1){
  h <- hist(pval3[ind3],main=paste(label," z=+1"),xlab="mixed predictive p-value",xlim=xlim,breaks=breaks)
  #mtext(sum(h$counts))
  #title(sub=paste("only genes prob(z)>",probz))
}
else plot(c(0,1),c(0,1),type="n")


########### qq-plots conditional on z, only genes where strong inference on zg

if(sum(ind1)>1){
  plot(sort(pval1[ind1]),ylab="",xlab=paste(sum(!is.nan(pval1[ind1]))," genes"),
       main="qq-plot, z=-1",ylim=c(0,1))  
  #abline(0,1/n1)
  abline( min.lhs , (1-min.lhs)/n1 )
  #title(sub=paste(sum(!is.nan(pval1[ind1]))," individuals"))
  #mtext(round(ks.test(pval1[ind1],"punif")$p.value,d=3))
}
else plot(c(0,1),c(0,1),type="n")

if(sum(ind2)>1){
  plot(sort(pval2[ind2]),ylab="",xlab=paste(sum(!is.nan(pval2[ind2]))," genes"),
       main="qq-plot, z=0",ylim=c(0,1))  
  abline(0,1/n2)
  #title(sub=paste(sum(!is.nan(pval2[ind2]))," individuals"))
  #mtext(round(ks.test(pval2[ind2],"punif")$p.value,d=3))
}
else plot(c(0,1),c(0,1),type="n")

if(sum(ind3)>1){
  plot(sort(pval3[ind3]),ylab="",xlab=paste(sum(!is.nan(pval3[ind3]))," genes"),
       main="qq-plot, z=+1",ylim=c(0,1))  
  #abline(0,1/n3)
  abline( 0 , max.rhs/n3 )
  #title(sub=paste(sum(!is.nan(pval3[ind3]))," individuals"))
  #mtext(round(ks.test(pval3[ind3],"punif")$p.value,d=3))
}
else plot(c(0,1),c(0,1),type="n")



}

