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
# fdr.cutoff -  which FDR cutoff(s) to use 

# OUTPUT 

# a data frame with the number of differentially expressed genes 
#              and the corresponding cutoffs on the tail posterior probability 

summaryTailPP <- function(tpp.res, fdr.cutoff =c(0.01, 0.05))
{
  ng = rep(NA, length(fdr.cutoff))
  pg = rep(NA, length(fdr.cutoff))

  for(k in 1:length(fdr.cutoff))
  {
   ng[k] <- sum(tpp.res$FDR <= fdr.cutoff[k], na.rm=T) 
   pg[k] <- min(tpp.res$tpp[tpp.res$FDR <= fdr.cutoff[k] & !is.na(tpp.res$FDR)])
  }
  output = data.frame(FDR.cutoff = paste(round(fdr.cutoff*100,1), "%", sep=""), Tailpp.cutoff=pg, Number.DEgenes = ng)

  return(output)
}
