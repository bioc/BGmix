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
# Calculate tail posterior probabality and FDR, plot a histogram

#TailPP <- function(res, nreps, params, paired=F, alpha=0.05, N = 5000, prec=0.05, p.cut = 0.7, plots  = T, plot.pi0=F)

# INPUT

#  samples from posterior distribution, genes - rows, iterations - columns:
# res - list object output from 'ccTrace'
# nreps - vector length 2 containing the number of replicates in each condition 
# params - list object output from 'ccParams'
# paired - logical. TRUE for paired design, FALSE for unpaired.
# alpha  - parameter of the tail posterior probability (1-alpha/2 quantile)
# N - simulation size for tail posterior probability under H0
# prec   - estimate CDF of tail posterior probability under H0 at points 1 - k*prec, k - integer
# p.cut  - to save time, calculate FDR only for cutoffs on tail posterior probability > p.cut 
# plots - if True, makes plots of the histogram of tail posterior probability with the null density  and of FDR
# plot.pi0  - if True, diagnostic plot the estimated pi0 at different locations and the median estimate 

# OUTPUT: 

# tpp - vector of tail posterior probabilities with parameter alpha, one per gene
# pi0 - estimate of pi0
# FDR - (smoothed) estimate of FDR for all (distinct) cutoffs > p.cut

TailPP <- function(res, nreps, params, paired=F, alpha=0.05, N = 5000, prec=0.05, p.cut = 0.7, plots  = T, plot.pi0=F)
{

  
  if(paired){
    dif <- res$beta[,,]
    var1 <- res$sig2[,,]
    var2 <- NULL
    n.rep1 <- nreps[1]
    n.rep2 <- NULL
    a1 <- params$maa[1]
    a2 <- NULL
    print("CHECK THIS FUNCTION")
  }
  else{
    dif <- res$beta[2,,]
    var1 <- res$sig2[1,,]
    var2 <- res$sig2[2,,]
    n.rep1 <- nreps[1]
    n.rep2 <- nreps[2]
    a1 <- params$maa[1]
    a2 <- params$maa[2]
  }
  
  tpp <- TailPPCalculate(dif, var1, var2, n.rep1, n.rep2, alpha)

  if(is.null(n.rep2)) A = a1 + (n.rep1-1)/2 else A = max(a1 + (n.rep1-1)/2,  a2 + (n.rep2-1)/2) # or min: to make FDR est. more conservative
  pp0 = SimulateFromNull(A, N)
  FDR = FDRforTailPP(tpp, a1, a2, n.rep1, n.rep2, prec, p.cut, N, pp0, plot = plot.pi0)
 
  tpp.res = list(tpp = tpp, FDR=FDR$FDR, pi0 = FDR$pi0, pp0=pp0)

  if(plots)
  {
   par(mfrow=c(2,2))
   histTailPP(tpp.res)
   FDRplotTailPP(tpp.res, plot.TP = T)
   FDRplotTailPP(tpp.res, plot.TP = F)
   par(mfrow=c(1,1))
  }

  return(tpp.res)

}
