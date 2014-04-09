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
# Calculate FDR for the tail posterior probabality

# ------------------------------------------------------------------

# INPUT 

# tpp - vector of tail posterior probabilities
# a1  - posterior mean of the shape parameter of the inverse gamma distribution - prior for the variance in condition 1 
# a2  - posterior mean of the shape parameter of the inverse gamma distribution - prior for the variance in condition 2
# n.rep1 - number of replicates in condition 1
# n.rep2 - number of replicates in condition 2
# prec   - estimate CDF of tail posterior probability under H0 at points 1 - k*prec, k - integer
# p.cut  - to save time, calculate FDR only for cutoffs on tail posterior probability > p.cut 
# N - simulation size for tail posterior probability under H0
# plot  - if True, plot the estimated pi0 at different locations and the median estimate 

# OUTPUT: 

# pi0 - estimate of pi0
# FDR - estimate of FDR for all (distinct) cutoffs > p.cut
# -------------------------------------------------------

FDRforTailPP <- function(tpp, a1, a2=NULL, n.rep1, n.rep2=NULL, prec=0.05, p.cut = 0.7, N = 10000, pp0=NULL, plot = T)
{

 if(is.null(pp0))
 {
  if(is.null(n.rep2)) A = a1 + (n.rep1-1)/2 else A = max(a1 + (n.rep1-1)/2,  a2 + (n.rep2-1)/2) # or min: to make FDR est. more conservative
  pp0 = SimulateFromNull(A, N)
 }

 F0 = calcF0tail(tpp, A, h=prec, p.cut, N=N, pp0=pp0, F0.only=T) 

 pi0.est = EstimatePi0(tpp, pp0, plot)

# estimate FDR

 n = length(tpp)

 FDR = rep(NA, n)

 ord1 = rev(order(tpp))
 
 FDR[ord1] = pi0.est*(1-F0[ord1])/seq(1/n,1,1/n)

# smooth
 
 h0=max(min(diff(sort(tpp))), 0.005)

 x=tpp[rev(order(tpp))]
 y=FDR[rev(order(tpp))]
 x=x[! is.na(y)]
 y=y[! is.na(y)]*(1:length(x))

 FDR.ksmooth = locpoly(x, y, bandwidth=2*h0)

 h = diff(FDR.ksmooth$x)[3]
 FDRs = rep(NA,n)
 for(k in 0:ceiling((1-p.cut)/h))
  FDRs[tpp <= 1 - k*h & tpp > 1 - (k+1)*h] = mean(FDR.ksmooth$y[FDR.ksmooth$x <= 1- k*h & FDR.ksmooth$x > 1- (k+1)*h])

 FDRs[rev(order(tpp))] =  FDRs[rev(order(tpp))]/(1:n)

 return(list(pi0=pi0.est, FDR=FDRs))
}

