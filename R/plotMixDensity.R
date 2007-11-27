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

"plotMixDensity" <-
function(res, predres, ybar, ss){

############################################################

par(mfrow=c(2,3))

diff.obs <- ybar[,2]-ybar[,1]
diff.pred <- predres$ybar.pred3[2,,] - predres$ybar.pred3[1,,]

#qqplot(diff.obs,diff.pred)

ind1 <- res$pc[,1]>0.6
ind2 <- res$pc[,2]>0.6
ind3 <- res$pc[,3]>0.6
print(sum(ind1))
print(sum(ind2))
print(sum(ind3))

hist(log(ss[,1]),freq=F,nc=20,main="Sum squares, cond 1")
lines(density(log(predres$ss.pred2[1,,])),col="red")

hist(log(ss[,2]),freq=F,nc=20,main="Sum squares, cond 2")
lines(density(log(predres$ss.pred2[2,,])),col="red")


hist(diff.obs,freq=F,nc=20,main="log fold change")
lines(density(diff.pred),col="red")

hist(diff.obs[ind1],freq=F,nc=20,main="log fold change, z=-1")
lines(density(diff.pred[ind1]),col="red")

hist(diff.obs[ind2],freq=F,nc=20,main="log fold change, z=0")
lines(density(diff.pred[ind2]),col="red")

hist(diff.obs[ind3],freq=F,nc=20,main="log fold change, z=+1")
lines(density(diff.pred[ind3]),col="red")




}

