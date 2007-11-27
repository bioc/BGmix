#  This file is part of BGmix, a fully Bayesian model for
#  differential expression.
#  Copyright 2007 Natalia Bochkina <N.Bochkina@ed.ac.uk>
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

#========================================================================================

# estimate pi0 - proportion of non-differentially expressed 

# tpp -  observed tail posterior probability
# pp0 - a vector of tail posterior probability under H0
# plot  - if True, plot the estimated pi0 at different locations and the median estimate 

EstimatePi0 <- function(tpp, pp0, plot = T)
{

lambda = seq(0.05,0.5, 0.01)

pi0 = rep(NA, length(lambda))

for(i in 1:length(lambda))
 pi0[i] = mean(tpp < lambda[i])/mean(pp0 < lambda[i])


pi0.est = min(median(pi0[lambda<= 0.3 & lambda >= 0.2]), 1)

if(plot) { 
plot(lambda, pi0, ylim=c(0,1))
abline(h=pi0.est) 
}

return(pi0.est)

}
