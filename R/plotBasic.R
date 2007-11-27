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

"plotBasic" <-
function(res, ybar, ss, q.mean=T, q.diff=T, q.sig=T, q.volcano=T){
  
ntau <- dim(ss)[2]

#print("WARNING: These plots designed for alpha,delta parametrization.")

print("These plots are designed for differential expression model.")

#print(dim(res$msig2))
#print(dim(ss))

##if(names(dev.cur())=="X11") X11()
par(mfrow=c(2,3))

## overall exp
if(q.mean) plotCompare(0.5*(ybar[,1]+ybar[,2]),res$mbeta[,1],xlab="Data mean",ylab="alpha",main="alpha")
  
## log fold change
if(q.diff) plotCompare(ybar[,2]-ybar[,1],res$mbeta[,2],xlab="Data diff",ylab="delta",main="delta")

## variances
#if(q.sig){
#  for(i in 1:ntau){
#    plotCompare(ss[,i],res$msig2[,i],log="xy",xlab=paste("Data variance, condition",i),
#                ylab=paste("sigma^2, condition",i))
#  }
#}

plotCompare(ss[,1],res$msig2[,1],log="xy",xlab=paste("Data variance, condition",1),
                ylab=paste("sigma^2, condition",1),main=paste("sigma^2, condition",1))
  
## posterior probs v. delta
if(q.volcano){
  plot(res$mbeta[,2],res$pc[,1],xlab="delta",ylab="Prob(under-expressed)",main="Prob(under-expressed)")
  plot(res$mbeta[,2],res$pc[,2],xlab="delta",ylab="Prob(not differentially expressed)",main="Prob(not differentially expressed)")
  plot(res$mbeta[,2],res$pc[,3],xlab="delta",ylab="Prob(over-expressed)",main="Prob(over-expressed)")
}
  
}

