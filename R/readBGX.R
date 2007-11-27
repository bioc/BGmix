#  This file is part of BGmix, a fully Bayesian model for
#  differential expression.
#  Copyright 2007 Ernest Turro <e.turro@imperial.ac.uk>
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

readBGX <- function(path) {

   muave <- matrix(scan(file.path(path,"muave")), ncol=2)

   sigma1 <- matrix(scan(file.path(path, "sigma2.1")), ncol=1024)
   sigma2 <- matrix(scan(file.path(path, "sigma2.2")), ncol=1024)

   sigma2ave <- matrix(nrow=nrow(sigma1), ncol=2)
   sigma2ave[,1] <- apply(sigma1,1,mean)
   sigma2ave[,2] <- apply(sigma2,1,mean)

return(list(ybar=muave,ss=sigma2ave))

}
