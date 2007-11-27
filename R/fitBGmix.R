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

"fitBGmix" <-
function(ngenes, nconds, neffects, ncomps, ntau, niter, nburn, nthin, jstar, 
	seed, zg.init, nlam, move.choice.bz, move.choice.cut, move.choice.aa, move.choice.lam,
	    move.choice.tau, move.choice.eta, like.choice, trace.out, trace.pred,
	aa.init, gg, hh, tau.var, eta.up.init, eta.down.init, 
	lambda.up.init, lambda.down.init, aa.eta, bb.eta, lam1, lam2, nu1, nu2, 
	nu0, df, beta.init, tau.init, bb.init, sig.aa, tau.eps, ydata,
	ybar, ss, nreps, nrepstot, xx, indtau, datafilename.ybar, xfilename, itfilename, rundir){
  

# free allocated memory on user interrupt/end of simulation
on.exit(.C("freeBGmixMemory", as.integer(ngenes), as.integer(neffects), PACKAGE = "BGmix"))


aa <- .C("BGmix",
         as.integer(ngenes),
         as.integer(nconds),
         as.integer(neffects),
         as.integer(ncomps),
         as.integer(ntau),
         as.integer(niter),
         as.integer(nburn),
         as.integer(nthin),
         as.integer(jstar),
         as.integer(seed),
         as.integer(zg.init),
         as.integer(nlam),
         as.integer(move.choice.bz),
         as.integer(move.choice.cut),
         as.integer(move.choice.aa),
         as.integer(move.choice.lam),
         as.integer(move.choice.tau),
         as.integer(move.choice.eta),
         as.integer(like.choice),
         as.integer(trace.out),
         as.integer(trace.pred),
         as.double(aa.init),
         as.double(gg),
         as.double(hh),
         as.double(tau.var),
         as.double(eta.up.init),
         as.double(eta.down.init),
         as.double(lambda.up.init),
         as.double(lambda.down.init),
         as.double(aa.eta),
         as.double(bb.eta),
         as.double(lam1),
         as.double(lam2),
         as.double(nu1),
         as.double(nu2),
         as.double(nu0),
         as.double(df),
         as.double(beta.init),
         as.double(tau.init),
         as.double(bb.init),
         as.double(sig.aa),
         as.double(tau.eps),         
         as.double(ydata),
         as.double(ybar),
         as.double(ss),
         as.integer(nreps),
         as.integer(nrepstot),
         as.double(xx),
         as.integer(indtau),
         as.character(datafilename.ybar),
         as.character(xfilename),
         as.character(itfilename),
         as.character(paste(character(length=nchar(rundir)+16), collapse="x")), ## output from C code
         as.character(rundir),
         PACKAGE="BGmix"
)

return(aa[[53]])


}

