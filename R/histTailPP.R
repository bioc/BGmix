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

#------------------------------------------------------------------------
# plot on the output of function "TailPP" 
#------------------------------------------------------------------------
#------------------------------------------------------------------------

#INPUT

# tpp.res - output from function "TailPP": list containing tail posterior probability "tpp" and "FDR"
# bw - bandwidth of the density of distribution of tail posterior probability under H0 (parameter to "density" function)
# parameters below are for the "hist" function for the tail posterior probability:
#  xlim - limits of the x axis
#  nc - number of classes of the histogram


histTailPP <- function(tpp.res, bw=0.05, xlim=c(0,1),nc=10)
{
 hs = hist(tpp.res$tpp, probability=T, xlim=c(0,1),plot=F, nc=nc)
 ymax=max(hs$density[hs$breaks[1:(length(hs$breaks)-1)] >= xlim[1]])
 hist(tpp.res$tpp, probability=T, xlim=xlim, ylim=c(0, ymax),xlab="Tail posterior probability", main = "Histogtam of tail pp", nc=nc)
 
# plot density of tail posterior probabilities under H0

 lines(density(tpp.res$pp0, bw=bw, from=0,to=1))

}
