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
# Calculate tail posterior probabality

# TailPP.Calculate <- function(diff, var1, var2 = NULL, n.rep1, n.rep2=NULL, alpha=0.05)

# INPUT

#  samples from posterior distribution, genes - rows, iterations - columns:
# diff - of the difference, 
# var1 - variance in condition 1
# var2 - variance in condition 2 (NULL if one condition) 
# n.rep1 - number of replicates in condition 1
# n.rep2 - number of replicates in condition 2
# alpha  - parameter of the tail posterior probability (1-alpha/2 quantile)

# OUTPUT

# vector of tail posterior probabilities with parameter alpha, one per gene
#----------------------------------------------------------------------

TailPPCalculate <- function(diff, var1, var2 = NULL, n.rep1, n.rep2=NULL, alpha=0.05)
{

# number of genes 
 n.genes = dim(diff)[1]

# number of MCMC Citerations  
 n.it = dim(diff)[2]

# variance of diff 
 
 if(is.null(n.rep2)) vart = var1/n.rep1 else vart = var1/n.rep1+var2/n.rep2

# t parameter

 tp = diff/sqrt(vart)

# threshold on t
 ta = qnorm(1-alpha/2) 

 tps = array(0,dim=dim(tp))
 tps[abs(tp) > ta]=1

 tpp = rowsum(t(tps), rep(1,n.it))/n.it
 
 return(tpp[1,])

}
