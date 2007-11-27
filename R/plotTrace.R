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

"plotTrace" <-
function(res, q.beta=T, q.sig=T, q.z=T, ind.genes=(1:3)){

summ <- res$summ

move.choice.aa <- summ$move.choice.aa
move.choice.lam <- summ$move.choice.lam
move.choice.eta <- summ$move.choice.eta
  
  
### options
q.aa <- move.choice.aa==1
if(is.null(move.choice.eta)) q.eta <- FALSE
else q.eta <- move.choice.eta==1
if(is.null(move.choice.lam)) q.lam <- FALSE
else q.lam <- move.choice.lam==1


################# scalars

par(mfrow=c(4,1))

if(q.aa){
  plot(res$aa[,1],type="l",main="a[1]")
  plot(res$aa[,2],type="l",main="a[2]")
}
plot(res$bb[,1],type="l",main="b[1]")
plot(res$bb[,2],type="l",main="b[2]")

if(names(dev.cur())=="X11") X11()
par(mfrow=c(4,1))

if(q.eta){
  plot(res$eta[1,],type="l",main="eta[1]")
  plot(res$eta[2,],type="l",main="eta[2]")
}
if(q.lam){
  plot(res$lam[1,],type="l",main="lam[1]")
  plot(res$lam[2,],type="l",main="lam[2]")
}

if(names(dev.cur())=="X11") X11()
par(mfrow=c(4,1))

plot(res$wtc[,1],type="l",main="w[1]")
plot(res$wtc[,2],type="l",main="w[2]")
plot(res$wtc[,3],type="l",main="w[3]")

################# large matrices

for(i in 1:length(ind.genes)){
  if(q.beta){
    plot(res$beta[1,ind.genes[i],],type="l")
    plot(res$beta[2,ind.genes[i],],type="l")
  }
  if(q.sig){
    plot(res$sig2[1,ind.genes[i],],type="l",log="y")
    plot(res$sig2[2,ind.genes[i],],type="l",log="y")
  }
  if(q.z){
    plot(res$zg[,ind.genes[i]],type="l")
  }
}
  
}

