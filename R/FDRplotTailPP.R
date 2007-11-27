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
# nmax  -   plot FDR/TP for the gene lists up to size "nmax"
# plot.TP - if true, plot Number of True Positives, otherwise plot FDR

FDRplotTailPP <- function(tpp.res, nmax = sum(! is.na(tpp.res$FDR)), plot.TP = F)
{
 ord = rev(order(tpp.res$tpp))
 n <- length(tpp.res$tpp)
 
 if(plot.TP) plot(tpp.res$FDR[ord]*(1:n), xlab="Number of genes in the list", ylab="Number of FP", xlim=c(0,nmax), type='l', main="Estimated number of FP") 
  else plot(tpp.res$FDR[ord], xlab="Number of genes in the list", ylab="FDR", xlim=c(0,nmax), type='l', main = "Estimated FDR")

}

