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

"ccSummary" <-
function(filedir){

file <- file.path(filedir,"summary.txt")
summ <- readLines(file)


i <- agrep("genes",summ,max=0)
ngenes <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(ngenes))

i <- agrep("conds",summ,max=0)
nconds <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(nconds))

i <- agrep("effects",summ,max=0)
neffects <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(neffects))

i <- agrep("mix comps",summ,max=0)
ncomps <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(ncomps))

i <- agrep("variance",summ,max=0)
ntau <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(ntau))

i <- agrep("mix prior",summ,max=0)
jstar <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(jstar))



i <- agrep("niter",summ,max=0)
niter <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(niter))

i <- agrep("nburn",summ,max=0)
nburn <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(nburn))

i <- agrep("nthin",summ,max=0)
nthin <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(nthin))

i <- agrep("seed",summ,max=0)
seed <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(seed))

i <- agrep("move",summ,max=0)
move.choice.bz <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(move.choice.bz))

i <- agrep("fullBayes",summ,max=0)
move.choice.cut <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(move.choice.cut))

i <- agrep("aa updated",summ,max=0)
move.choice.aa <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(move.choice.aa))

i <- agrep("eta updated",summ,max=0)
#print(i)
#print(length(strsplit(summ[i]," ")))
if(length(strsplit(summ[i]," "))>0) move.choice.eta <- as.numeric(strsplit(summ[i]," ")[[1]][1])
else move.choice.eta <- NULL
#print(as.numeric(move.choice.eta))

i <- agrep("lambdas updated",summ,max=0)
move.choice.lam <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(move.choice.lam))

i <- agrep("tau's",summ,max=0)
#print(i)
#print(length(strsplit(summ[i]," ")))
if(length(strsplit(summ[i]," "))>0) move.choice.tau <- as.numeric(strsplit(summ[i]," ")[[1]][1])
else move.choice.tau <- NULL
#print(as.numeric(move.choice.tau))

i <- agrep("t-dist",summ,max=0)
like.choice <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(like.choice))

i <- agrep("(params)",summ,max=0)
trace.out <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(trace.out))

i <- agrep("(predictive)",summ,max=0)
trace.pred <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(trace.pred))


i <- agrep("eta_up",summ,max=0)
eta.up.init <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(eta.up.init))

i <- agrep("eta_down",summ,max=0)
eta.down.init <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(eta.down.init))

i <- agrep("lambda+ for",summ,max=0)
lambda.up.init <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(lambda.up.init))

i <- agrep("lambda- for",summ,max=0)
lambda.down.init <- as.numeric(strsplit(summ[i]," ")[[1]][1])
#print(as.numeric(lambda.down.init))




#print(strsplit(summ[i]," "))



return(list(ngenes=ngenes,nconds=nconds,neffects=neffects,ncomps=ncomps,ntau=ntau,jstar=jstar,
       niter=niter,nburn=nburn,nthin=nthin,seed=seed,
       move.choice.bz=move.choice.bz,move.choice.cut=move.choice.cut,move.choice.aa=move.choice.aa,
       move.choice.eta=move.choice.eta,move.choice.lam=move.choice.lam,move.choice.tau=move.choice.tau,
       like.choice=like.choice,trace.out=trace.out,trace.pred=trace.pred,
            lambda.up.init=lambda.up.init,lambda.down.init=lambda.down.init,
            eta.up.init=eta.up.init,eta.down.init=eta.down.init))

}

