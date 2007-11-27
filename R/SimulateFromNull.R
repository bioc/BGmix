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


SimulateFromNull <- function(A, N=10000, N1 = 1000)
{
 x0 = rt(N, 2*A)
 
 x01 = matrix(x0, ncol=N1, nrow=N)

 zz = matrix(sqrt(rgamma(N*N1,A,1)/A), ncol=N1, nrow=N)

 xx=t(x01*zz)
 pp0 = as.vector(rowsum(pnorm(-2-xx) + pnorm(-2+xx), rep(1,N1)))/N1

 return(pp0)
}
