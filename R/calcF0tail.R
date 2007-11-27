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

#=================================================
# Auxiliary functions 
#=================================================

# -------------------------------------------
# input: pp- tail posterior  probability

calcF0tail <- function(pp, A,  h=0.02, p.cut=0.6, N=100000, pp0=SimulateFromNull(A, N), correct.zero = F, F0.only=T)
{

 if(is.null(dim(pp))) pp = as.matrix(pp)
 n = dim(pp)[1]
 K = dim(pp)[2]

 f0 = array(dim=c(n,K))
 F0 = array(dim=c(n,K))

 for(k in 1:K) 
 {

  seq1 = (1:n)[pp[,k] > p.cut]
  for(i in seq1)
  {
   if(! F0.only) f0[i,k] = mean(pp0 >= pp[i,k]-h/2 & pp0 <= pp[i,k]+h/2)/h
   F0[i,k] = mean(pp0 < pp[i,k])
  } 

  
  if(correct.zero & ! F0.only)
  {
   f0.summ = tapply(f0[,k], INDEX = pp[,k], mean)
   pp.summ = tapply(pp[,k], INDEX = pp[,k], mean)
#   plot(pp.summ, f0.summ)

   seq01 = match(pp[pp[,k] < 0.05,k], pp.summ)

   f0[pp < 0.05,k] = f0.summ[12 - seq01]
  }
 }

 if(F0.only) return(F0) else return(list(f0=f0, F0=F0))

}
