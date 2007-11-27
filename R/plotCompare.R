#  This file is part of BGmix, a fully Bayesian model for
#  differential expression.
#  Copyright 2007 Alex Lewin <a.m.lewin@imperial.ac.uk>
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

"plotCompare" <-
function(var1,var2,limi=0,
xlab=substitute(var1),ylab=substitute(var2),
log="",title="",...){

if(length(limi)==1) limi <- c(min(var1,var2,na.rm=T),max(var1,var2,na.rm=T))
plot.default(var1,var2,xlim=limi,ylim=limi,xlab=xlab,ylab=ylab,log=log,...)
title(title)
abline(0,1)

#print("")
#print(paste(xlab,ylab))
#print("bias")
#print( mean(var2-var1) )
#print("rms")
#print( sqrt(mean((var2-var1)**2)) )

return(limi)

}

