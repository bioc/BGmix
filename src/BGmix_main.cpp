/*
*  This file is part of BGmix, a fully Bayesian model for
*  differential expression.
*  Copyright 2007 Alex Lewin <a.m.lewin@imperial.ac.uk>
*  Thanks to Ernest Turro for the R stuff
*
*  BGmix is free software; you can redistribute it and/or modify it
*  under the terms of the GNU General Public License, version 2, as
*  published by the Free Software Foundation.
*
*  BGmix is distributed in the hope that it will be useful, but
*  WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*  General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include "updates_BGmix.h"
#include "rand.hh"
#include "rundir.hh"
#include "BGmix_main.h"

#ifdef USING_R
 extern "C" {
  #include <R.h> // for flushing console, allowing user interrupts, and printing to console
  #if ( defined(HAVE_AQUA) || defined(WIN32) )
    #define FLUSH {R_FlushConsole(); R_ProcessEvents();}
  #else
    #define FLUSH R_FlushConsole();
  #endif
  #define PRINTF(...) Rprintf((char*) __VA_ARGS__)
/* Register routines, allocate resources. */
 #define R_NO_REMAP
  #include <Rinternals.h>
  #include <R_ext/Rdynload.h>
 }
static const R_CMethodDef cMethods[] = {
    {"BGmix", (DL_FUNC)&BGmix, 54},
    {"freeBGmixMemory", (DL_FUNC)&freeBGmixMemory, 2},
    {NULL, NULL, 0}
  };
  void R_init_BGmix(DllInfo *info) {
    R_registerRoutines(info, cMethods,NULL,NULL,NULL);
  }
  void R_unload_BGmix(DllInfo *info) { } 
#else
  #define FLUSH fflush(stdout);
  #define PRINTF(...) printf(__VA_ARGS__)
  #define CARRIAGERETURN "\r"
#endif

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

using namespace std;

void stringcpy(char* p, const std::string& s)
{
  s.copy(p,std::string::npos);
  p[s.length()] = 0; // Add C-style string terminator
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

// The following are global so they can be deletd in freeBGmixMemory 
// if there is a user interrupt in R

static double **ybar, **ss, **ydata, **beta, **tau, **gamma1,
  **mean_beta, **mean_tau, **meanlogscale_tau, **mean_sig2, **mean_gamma1,
  **meanlogscale_gamma1, *mean_zg, *prob0, *prob1, *prob2, 
  **ybar_pred1, **ybar_pred2, **ybar_pred3, **ybar_pred4, **ss_pred1, **ss_pred2,
  **pval_post_ss, **pval_mix_ss, **pval_partial_ss, 
  **pval_post_ybar, **pval_mix1_ybar, **pval_mix2_ybar, **pval_mix3_ybar, *pval_partial_ybar, 
  **norm_ss, *norm_ybar;
static double *aa, *bb, *nu, *wtc, *df, **xx, *mean_aa, *mean_bb, *mean_wtc;
static int *nalloc;
static int *zg;

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

void BGmix(int &ngenes, int &nconds, int &neffects, int &ncomps, int &ntau, int &niter, 
	    int &nburn, int &nthin, int &jstar, int &seed, int &zg_init, int &nlam, 
	    int &move_choice_bz, int &move_choice_cut, int &move_choice_aa, int &move_choice_lam,
	    int &move_choice_tau, int &move_choice_eta, int &like_choice, int &trace_out, int &trace_pred,
	    double &aa_const, double &gg_const, double &hh_const, double &tau_var_const,
	    double &eta_up_const, double &eta_down_const, double &lambda_up_const, double &lambda_down_const, 
	    double &aa_eta_const, double &bb_eta_const, double &lam1, double &lam2,
	    double &nu1_const, double &nu2_const, double &nu0_const, double &df_const,
	    double &beta_init, double &tau_init, double &bb_init, double &sig_aa, double &tau_eps,
	    double* ydata_in, double* ybar_in, double* ss_in, int* nreps, int &nrepstot, 
	   double* xx_in, int* indtau, char** dataname, char** xname, char** indtauname,
	   char** dirname, char **basepath)
{

//////////////////////////// make vectors into matrices ///////////////////////////////////

ybar = new double*[ngenes];
for(int i=0; i<ngenes; ++i) ybar[i] = new double[nconds];
ss = new double*[ngenes];
for(int i=0; i<ngenes; ++i) ss[i] = new double[nconds];

ydata = new double*[ngenes];
for(int i=0; i<ngenes; ++i) ydata[i] = new double[nrepstot];

xx = new double*[neffects];
for(int i=0; i<neffects; ++i) xx[i] = new double[nconds];

for(int i=0; i<ngenes; ++i) {
  for(int j=0; j<nconds; ++j) {
	ybar[i][j] = ybar_in[i+ngenes*j];
	ss[i][j] = ss_in[i+ngenes*j];
	//if(i==2) cerr << ybar[i][j] << " ybar[i=2] " << endl;
	//if(i==2) cerr << ss[i][j] << " ss[i=2] " << endl;
  }
  for(int j=0; j<nrepstot; ++j) {
	ydata[i][j] = ydata_in[i+ngenes*j];
	//if(i==2) cerr << ydata[i][j] << " ydata[i=2] " << endl;
  }
}

for(int i=0; i<neffects; ++i) {
  for(int j=0; j<nconds; ++j) {
	xx[i][j] = xx_in[i+neffects*j];
	//cerr << xx[i][j] << " xx[i][j] " << endl;
}}

//for(int j=0; j<nconds; ++j) {
//  cerr << nreps[j] << " nreps[j]" << endl;
//  cerr << indtau[j] << " indtau[j]" << endl;
//}

// cerr << "nrepstot = " << nrepstot << endl;

///////////////////////////////////////////// init vals ///////////////////////////////////

beta = new double*[ngenes];
for(int i=0; i<ngenes; ++i) beta[i] = new double[neffects];
tau = new double*[ngenes];
for(int i=0; i<ngenes; ++i) tau[i] = new double[ntau];
gamma1 = new double*[ngenes];
for(int i=0; i<ngenes; ++i) gamma1[i] = new double[nrepstot];

aa = new double[ntau];
bb = new double[ntau];
double gg,hh,tau_var;
zg = new int[ngenes];
double eta_down,eta_up,lambda_up,lambda_down,aa_eta,bb_eta;
nu = new double[ncomps];
wtc = new double[ncomps];
nalloc = new int[ncomps];
df = new double[nconds];

//// init stochastic params
 for(int i=0; i<ngenes; ++i) {
   for(int j=0; j<neffects; ++j) {
     if(neffects==2 & nconds==2){
       beta[i][0] = xx[1][1]*ybar[i][0] - xx[1][0]*ybar[i][1];
       beta[i][1] = xx[0][0]*ybar[i][1] - xx[0][1]*ybar[i][0];
     }
     else beta[i][j] = beta_init;
     //beta[i][j] = beta_init;
   }
   for(int j=0; j<ntau; ++j) {
     tau[i][j] = tau_init;
   }
   for(int j=0; j<nrepstot; ++j) {
     gamma1[i][j] = 1;
   }
   if(jstar!=-1 & ncomps==3){
     if(beta[i][jstar]<0.01 & beta[i][jstar]>-0.01) zg[i]=1;
     else if(beta[i][jstar]>0) zg[i]=2;
     else zg[i]=0;
   }
   else zg[i] = zg_init;
   //zg[i] = zg_init;
 }
 
for(int j=0; j<ntau; ++j) bb[j] = bb_init;

//// init nalloc and mix component weights
for(int j=0; j<ncomps; ++j) nalloc[j] = 0;
for(int i=0; i<ngenes; ++i) {
	if(zg[i]==0) nalloc[0]++;
	else if(zg[i]==2)	nalloc[2]++;
	else nalloc[1]++;
}
for(int j=0; j<ncomps; ++j) wtc[j] = nalloc[j]/float(ngenes);

// constants for variances here 
for(int j=0; j<ntau; ++j) aa[j] = aa_const;
gg = gg_const;
hh = hh_const;
tau_var = tau_var_const;

// df for t
for(int j=0; j<nconds; ++j) df[j] = df_const;

// other params used in mixture components (NB here ncomps=3 always)
// eta_down,eta_up sometimes fixed,sometimes stoch
eta_down = eta_down_const;
eta_up = eta_up_const;
lambda_up = lambda_up_const;
lambda_down = lambda_down_const;
aa_eta = aa_eta_const;
bb_eta = bb_eta_const;
nu[0] = nu0_const;
nu[1] = nu1_const;
nu[2] = nu2_const;


///////////////////////////////////////////// init accumulators ///////////////////////////////////

mean_beta = new double*[ngenes];
for(int i=0; i<ngenes; ++i) mean_beta[i] = new double[neffects];

mean_zg = new double[ngenes];
prob0 = new double[ngenes];
prob1 = new double[ngenes];
prob2 = new double[ngenes];

for(int i=0; i<ngenes; ++i) {
 mean_zg[i] = 0;
 prob0[i] = 0;
 prob1[i] = 0;
 prob2[i] = 0;
 for(int j=0; j<neffects; ++j) {
	mean_beta[i][j] = 0;
}}

mean_sig2 = new double*[ngenes];
for(int i=0; i<ngenes; ++i) mean_sig2[i] = new double[ntau];
mean_tau = new double*[ngenes];
for(int i=0; i<ngenes; ++i) mean_tau[i] = new double[ntau];
meanlogscale_tau = new double*[ngenes];
for(int i=0; i<ngenes; ++i) meanlogscale_tau[i] = new double[ntau];

 for(int i=0; i<ngenes; ++i) {
   for(int j=0; j<ntau; ++j) {
     mean_sig2[i][j] = 0;
     mean_tau[i][j] = 0;
     meanlogscale_tau[i][j] = 0;
   }
 }

  mean_gamma1 = new double*[ngenes];
  for(int i=0; i<ngenes; ++i) mean_gamma1[i] = new double[nrepstot];
  meanlogscale_gamma1 = new double*[ngenes];
  for(int i=0; i<ngenes; ++i) meanlogscale_gamma1[i] = new double[nrepstot];

  for(int i=0; i<ngenes; ++i) {
    for(int j=0; j<nrepstot; ++j) {
      mean_gamma1[i][j] = 0;
      meanlogscale_gamma1[i][j] = 0;
    }
  }

// don't need init cos not part of MCMC chain
ybar_pred1 = new double*[ngenes];
for(int i=0; i<ngenes; ++i) ybar_pred1[i] = new double[nconds];
ybar_pred2 = new double*[ngenes];
for(int i=0; i<ngenes; ++i) ybar_pred2[i] = new double[nconds];
ybar_pred3 = new double*[ngenes];
for(int i=0; i<ngenes; ++i) ybar_pred3[i] = new double[nconds];
ybar_pred4 = new double*[ngenes];
for(int i=0; i<ngenes; ++i) ybar_pred4[i] = new double[nconds];
ss_pred1 = new double*[ngenes];
for(int i=0; i<ngenes; ++i) ss_pred1[i] = new double[nconds];
ss_pred2 = new double*[ngenes];
for(int i=0; i<ngenes; ++i) ss_pred2[i] = new double[nconds];


pval_post_ss = new double*[ngenes];
for(int i=0; i<ngenes; ++i) pval_post_ss[i] = new double[nconds];
pval_mix_ss = new double*[ngenes];
for(int i=0; i<ngenes; ++i) pval_mix_ss[i] = new double[nconds];
pval_partial_ss = new double*[ngenes];
for(int i=0; i<ngenes; ++i) pval_partial_ss[i] = new double[nconds];

/// now predicting delta but conditioned on zg
pval_post_ybar = new double*[ngenes];
for(int i=0; i<ngenes; ++i) pval_post_ybar[i] = new double[ncomps];
pval_mix1_ybar = new double*[ngenes];
for(int i=0; i<ngenes; ++i) pval_mix1_ybar[i] = new double[ncomps];
pval_mix2_ybar = new double*[ngenes];
for(int i=0; i<ngenes; ++i) pval_mix2_ybar[i] = new double[ncomps];
pval_mix3_ybar = new double*[ngenes];
for(int i=0; i<ngenes; ++i) pval_mix3_ybar[i] = new double[ncomps];
pval_partial_ybar = new double[ngenes];

// norm for partial pvals
norm_ss = new double*[ngenes];
for(int i=0; i<ngenes; ++i) norm_ss[i] = new double[nconds];
norm_ybar = new double[ngenes];


for(int i=0; i<ngenes; ++i) {
for(int j=0; j<nconds; ++j) {
	pval_post_ss[i][j] = 0;
	pval_mix_ss[i][j] = 0;
	pval_partial_ss[i][j] = 0;
	norm_ss[i][j] = 0;
}
}
for(int i=0; i<ngenes; ++i) {
for(int j=0; j<ncomps; ++j) {
	pval_post_ybar[i][j] = 0;
	pval_mix1_ybar[i][j] = 0;
	pval_mix2_ybar[i][j] = 0;
	pval_mix3_ybar[i][j] = 0;
}
pval_partial_ybar[i] = 0;
norm_ybar[i] = 0;
}

mean_bb = new double[ntau];
mean_aa = new double[ntau];
mean_wtc = new double[ncomps];

for(int j=0; j<ntau; ++j) mean_bb[j] = 0;
for(int j=0; j<ntau; ++j) mean_aa[j] = 0;
for(int j=0; j<ncomps; ++j) mean_wtc[j] = 0;

int n_acc=0;
int n_try=0;
int n_aa_acc=0;
int n_aa_try=0;
int n_tau_lognorm_acc=0;
int n_tau_lognorm_try=0;

double mean_eta_up=0,mean_eta_down=0;
double mean_lambda_up=0,mean_lambda_down=0;

double deviance1=0,deviance2=0,mean_dev1=0,mean_dev2=0;

////////////////////////////////////////// open files ///////////////////////////////////

std::string filename;

//std::string run_dir = rundir("run");
//stringcpy(*dirname,run_dir);

  // Create output file and directory
  char *tmpstr = new char[strlen(*basepath)+4];
  strcpy(tmpstr, *basepath); // Prepend basepath
  strcat(tmpstr,"/run");
  std::string run_dir = rundir(tmpstr);
  delete [] tmpstr;
  // In R, dirname is set to the ouput directory
  stringcpy(*dirname,run_dir);


filename=run_dir+"/summary.txt";
std::ofstream summary_file(filename.c_str());

filename=run_dir+"/trace_beta.txt";
std::ofstream beta_file(filename.c_str());
filename=run_dir+"/trace_sig2.txt";
std::ofstream sig2_file(filename.c_str());
filename=run_dir+"/trace_bb.txt";
std::ofstream bb_file(filename.c_str());
filename=run_dir+"/trace_aa.txt";
std::ofstream aa_file(filename.c_str());
filename=run_dir+"/trace_wtc.txt";
std::ofstream wtc_file(filename.c_str());
filename=run_dir+"/trace_zg.txt";
std::ofstream zg_file(filename.c_str());
filename=run_dir+"/trace_eta.txt";
std::ofstream eta_file(filename.c_str());
filename=run_dir+"/trace_lambda.txt";
std::ofstream lambda_file(filename.c_str());

filename=run_dir+"/trace_ybar_pred1.txt";
std::ofstream ybar_pred1_file(filename.c_str());
filename=run_dir+"/trace_ybar_pred2.txt";
std::ofstream ybar_pred2_file(filename.c_str());
filename=run_dir+"/trace_ybar_pred3.txt";
std::ofstream ybar_pred3_file(filename.c_str());
filename=run_dir+"/trace_ybar_pred4.txt";
std::ofstream ybar_pred4_file(filename.c_str());
filename=run_dir+"/trace_ss_pred1.txt";
std::ofstream ss_pred1_file(filename.c_str());
filename=run_dir+"/trace_ss_pred2.txt";
std::ofstream ss_pred2_file(filename.c_str());

filename=run_dir+"/mean_beta.txt";
std::ofstream mbeta_file(filename.c_str());
filename=run_dir+"/mean_sig2.txt";
std::ofstream msig2_file(filename.c_str());
filename=run_dir+"/mean_bb.txt";
std::ofstream mbb_file(filename.c_str());
filename=run_dir+"/mean_aa.txt";
std::ofstream maa_file(filename.c_str());
filename=run_dir+"/mean_tau.txt";
std::ofstream mtau_file(filename.c_str());
filename=run_dir+"/mean_wtc.txt";
std::ofstream mwtc_file(filename.c_str());
filename=run_dir+"/mean_zg.txt";
std::ofstream mzg_file(filename.c_str());
filename=run_dir+"/mean_eta.txt";
std::ofstream m_eta_file(filename.c_str());
filename=run_dir+"/mean_lambda.txt";
std::ofstream m_lambda_file(filename.c_str());

filename=run_dir+"/prob_class.txt";
std::ofstream probclass_file(filename.c_str());
filename=run_dir+"/pval_post_ss.txt";
std::ofstream post_ss_file(filename.c_str());
filename=run_dir+"/pval_mix_ss.txt";
std::ofstream mix_ss_file(filename.c_str());
filename=run_dir+"/pval_partial_ss.txt";
std::ofstream part_ss_file(filename.c_str());
filename=run_dir+"/pval_post_ybar.txt";
std::ofstream post_ybar_file(filename.c_str());
filename=run_dir+"/pval_mix1_ybar.txt";
std::ofstream mix1_ybar_file(filename.c_str());
filename=run_dir+"/pval_mix2_ybar.txt";
std::ofstream mix2_ybar_file(filename.c_str());
filename=run_dir+"/pval_mix3_ybar.txt";
std::ofstream mix3_ybar_file(filename.c_str());
filename=run_dir+"/pval_partial_ybar.txt";
std::ofstream part_ybar_file(filename.c_str());


////////////////////////////////// write summary file ///////////////////////////////////

summary_file << run_dir << " Output Directory" << std::endl;

summary_file << ngenes << " genes" << std::endl;
summary_file << nconds << " conds" << std::endl;
summary_file << neffects << " effects" << std::endl;
summary_file << ncomps << " mix comps" << std::endl;
summary_file << ntau << " variances" << std::endl;
summary_file << jstar << " = effect with mix prior (starts from 0)" << std::endl << std::endl;

summary_file << *dataname << " data file (ybar)" << std::endl;
summary_file << *xname << " x file" << std::endl;
summary_file << *indtauname << " indtau file" << std::endl << std::endl;

summary_file << niter << " niter" << std::endl;
summary_file << nburn << " nburn" << std::endl;
summary_file << nthin << " nthin" << std::endl;
summary_file << seed << " random seed" << std::endl;
summary_file << move_choice_bz << " beta,z prior/move choice" << std::endl;
summary_file << move_choice_cut << " 1=fullBayes, else cut" << std::endl;
summary_file << move_choice_aa << " 1=aa updated, else not" << std::endl;
summary_file << move_choice_eta << " 1=eta updated, else not" << std::endl;
summary_file << move_choice_lam << " 1=lambdas updated, else not" << std::endl;
summary_file << move_choice_tau << " 1=tau's are Gamma, 2= logNorm" << std::endl;
summary_file << like_choice << " 1=Normal,2=t-distn" << std::endl;
summary_file << trace_out << " 1=output trace (params), else not" << std::endl;
summary_file << trace_pred << " 1=output trace (predictive), else not" << std::endl << std::endl;

summary_file << aa_const << " a[c]" << std::endl;
summary_file << gg_const << " g (fixed)" << std::endl;
summary_file << hh_const << " h (fixed)" << std::endl;
summary_file << tau_var_const << " tau_var (fixed) for logNormal tau_gs" << std::endl;
summary_file << eta_up_const << " eta_up (fixed or init val)" << std::endl;
summary_file << eta_down_const << " eta_down (fixed or init val)" << std::endl;
summary_file << lambda_up_const << " lambda+ for Gam priors (fixed or init val)" << std::endl;
summary_file << lambda_down_const << " lambda- for Gam priors (fixed or init val)" << std::endl;
summary_file << aa_eta_const << " a_eta for Gam priors (fixed)" << std::endl;
summary_file << bb_eta_const << " b_eta for Gam priors (fixed)" << std::endl;
summary_file << lam1 << ", " << lam2 << ", " << nlam << " lam1,lam2,nlam for Gam priors (fixed)" << std::endl;
summary_file << nu0_const << " nu0 (fixed)" << std::endl;
summary_file << nu1_const << " nu1 (fixed)" << std::endl;
summary_file << nu2_const << " nu2 (fixed)" << std::endl;
summary_file << df_const << " df (fixed)" << std::endl;
summary_file << tau_eps << " tau_epsilon (fixed)" << std::endl << std::endl;

summary_file << beta_init << " initial beta" << std::endl;
summary_file << tau_init << " initial tau" << std::endl;
summary_file << bb_init << " initial bb" << std::endl;
summary_file << zg_init << " initial zg" << std::endl;
summary_file << sig_aa << " sig for RW update of aa" << std::endl;



///////////////////////////////////////////// updates ///////////////////////////////////

// Random is a class of random generators of various distributions (use eg. rand.Uniform() )
Random rand(seed);

//std::cerr << std::endl << "Burn-in: " << std::endl;
 PRINTF(" Burn-in: \n");

for(int iter=1; iter<=nburn; ++iter) {

  if(iter%(nburn/100)==0){
    //std::cerr << iter << " "; 
    PRINTF("%d ",iter);
    FLUSH;
  }

#ifdef USING_R
	R_CheckUserInterrupt();
#endif
	
	update_beta0(beta, tau, gamma1, xx, indtau, ybar, ydata, like_choice,
            ngenes, nconds, nreps, neffects, jstar, rand, summary_file);

	/// tau are Gamma
	if(move_choice_tau==1){
	  if(move_choice_cut==1) update_tau(beta, tau, gamma1, xx, indtau, ybar, ss, ydata,
				      aa, bb, like_choice, ngenes, nconds, ntau, nreps, neffects, rand);
	  else update_tau_cut(tau, ss, indtau, aa, bb, ngenes, 
			      nconds, ntau, nreps, rand);
	  update_bb(tau, aa, bb, gg, hh, ngenes, nconds, ntau, rand);
	  if(move_choice_aa==1) update_aa(sig_aa, tau, aa, bb, gg, hh, 
					  n_aa_acc, n_aa_try, ngenes, nconds, ntau, rand);
	}
	/// tau are logNorm
	else{
	  update_tau_logNorm(beta, tau, gamma1, xx, indtau, ybar, ss, ydata,
			     aa, bb, n_tau_lognorm_acc, n_tau_lognorm_try, like_choice,
			     ngenes, nconds, ntau, nreps, neffects, rand);
	  update_bb_logNorm(tau, aa, bb, gg, hh, ngenes, nconds, ntau, rand);
	  update_aa_logNorm(sig_aa, tau, aa, bb, tau_var, ngenes, nconds, ntau, rand);
	}

	if(like_choice==2) update_gamma(beta, tau, gamma1, xx, indtau, ydata, df, 
            ngenes, nconds, ntau, nreps, neffects, rand);

	/// mixture bit
	if(jstar!=-1) {
	  if(move_choice_bz==1){ 
	    update_z_beta1_joint1(zg, wtc, nalloc, 
				  eta_up, eta_down, beta, tau, gamma1, xx, indtau, ybar, ydata, like_choice, 
				  ngenes, nconds, nreps, ncomps, neffects, jstar, rand);
	    if(move_choice_eta==1) update_eta_unif(eta_up, eta_down, aa_eta, bb_eta, nalloc, rand);
	  }
	  else if(move_choice_bz==2){
	    update_z_beta1_joint2(zg, beta, nalloc, n_acc, n_try, 
				  wtc, eta_up, eta_down, tau, gamma1, xx, indtau, ybar, ydata, like_choice,
				  ngenes, nconds, nreps, ncomps, neffects, jstar, rand);
	    if(move_choice_eta==1) update_eta_unif(eta_up, eta_down, aa_eta, bb_eta, nalloc, rand);
	  }
	  else if(move_choice_bz==3 | move_choice_bz==4){
	    update_z_beta1_joint3(zg, beta, nalloc, n_acc, n_try, wtc, eta_up, eta_down, 
				  lambda_up, lambda_down, tau, gamma1, xx, indtau, ybar, ydata, like_choice, 
				  ngenes, nconds, nreps, ncomps, neffects, jstar, rand);
	    if(move_choice_eta==1) update_eta(eta_up, eta_down, lambda_up, lambda_down, aa_eta, bb_eta, 
					      zg, beta, nalloc, ngenes, nconds, jstar, rand);
	    if(move_choice_lam==1) update_lambda(lambda_up, lambda_down, eta_up, eta_down, 
						 lam1, lam2, nlam, zg, beta, nalloc, ngenes, jstar, rand);
	  }
	  else  if(move_choice_bz==5){
	    update_z_beta1_joint4(zg, beta, nalloc, n_acc, n_try, wtc, tau_eps, eta_up, eta_down, 
				  lambda_up, lambda_down, tau, gamma1, xx, indtau, ybar, ydata, like_choice, 
				  ngenes, nconds, nreps, ncomps, neffects, jstar, rand);
	    if(move_choice_eta==1) update_eta(eta_up, eta_down, lambda_up, lambda_down, aa_eta, bb_eta, 
					      zg, beta, nalloc, ngenes, nconds, jstar, rand);
	    if(move_choice_lam==1) update_lambda(lambda_up, lambda_down, eta_up, eta_down, 
						 lam1, lam2, nlam, zg, beta, nalloc, ngenes, jstar, rand);
	  }	  
	  update_wtc(wtc, nalloc, nu, ncomps, rand);
	}
}

//std::cerr << std::endl << "Main up-dates: " << std::endl;
 PRINTF("\n Main up-dates: \n");

for(int iter=1; iter<=niter; ++iter) {

  if(iter%(niter/100)==0){
    //std::cerr << iter << " "; 
    PRINTF("%d ",iter);
    FLUSH;
  }

#ifdef USING_R
	R_CheckUserInterrupt();
#endif

	update_beta0(beta, tau, gamma1, xx, indtau, ybar, ydata, like_choice, 
            ngenes, nconds, nreps, neffects, jstar, rand, summary_file);

	/// tau are Gamma
	if(move_choice_tau==1){
	  if(move_choice_cut==1) update_tau(beta, tau, gamma1, xx, indtau, ybar, ss, ydata, 
				    aa, bb, like_choice, ngenes, nconds, ntau, nreps, neffects, rand);
	  else update_tau_cut(tau, ss, indtau, aa, bb, ngenes, nconds, ntau, nreps, rand);
	  update_bb(tau, aa, bb, gg, hh, ngenes, nconds, ntau, rand);
	  if(move_choice_aa==1) update_aa(sig_aa, tau, aa, bb, gg, hh, 
					  n_aa_acc, n_aa_try, ngenes, nconds, ntau, rand);
	}
	/// tau are logNorm
	else{
	  update_tau_logNorm(beta, tau, gamma1, xx, indtau, ybar, ss, ydata,
			     aa, bb, n_tau_lognorm_acc, n_tau_lognorm_try, like_choice,
			     ngenes, nconds, ntau, nreps, neffects, rand);
	  update_bb_logNorm(tau, aa, bb, gg, hh, ngenes, nconds, ntau, rand);
	  update_aa_logNorm(sig_aa, tau, aa, bb, tau_var, ngenes, nconds, ntau, rand);
	}


	if(like_choice==2) update_gamma(beta, tau, gamma1, xx, indtau, ydata, df, 
            ngenes, nconds, ntau, nreps, neffects, rand);

	if(jstar!=-1) {
	  if(move_choice_bz==1){ 
	    update_z_beta1_joint1(zg, wtc, nalloc, 
				  eta_up, eta_down, beta, tau, gamma1, xx, indtau, ybar, ydata, like_choice, 
				  ngenes, nconds, nreps, ncomps, neffects, jstar, rand);
	    if(move_choice_eta==1) update_eta_unif(eta_up, eta_down, aa_eta, bb_eta, nalloc, rand);
	  }
	  else if(move_choice_bz==2){
	    update_z_beta1_joint2(zg, beta, nalloc, n_acc, n_try, 
				  wtc, eta_up, eta_down, tau, gamma1, xx, indtau, ybar, ydata, like_choice,
				  ngenes, nconds, nreps, ncomps, neffects, jstar, rand);
	    if(move_choice_eta==1) update_eta_unif(eta_up, eta_down, aa_eta, bb_eta, nalloc, rand);
	  }
	  else if(move_choice_bz==3 | move_choice_bz==4){
	    update_z_beta1_joint3(zg, beta, nalloc, n_acc, n_try, wtc, eta_up, eta_down, 
				  lambda_up, lambda_down, tau, gamma1, xx, indtau, ybar, ydata, like_choice, 
				  ngenes, nconds, nreps, ncomps, neffects, jstar, rand);
	    if(move_choice_eta==1) update_eta(eta_up, eta_down, lambda_up, lambda_down, aa_eta, bb_eta, 
					      zg, beta, nalloc, ngenes, nconds, jstar, rand);
	    if(move_choice_lam==1) update_lambda(lambda_up, lambda_down, eta_up, eta_down, 
						 lam1, lam2, nlam, zg, beta, nalloc, ngenes, jstar, rand);
	  }
	  else  if(move_choice_bz==5){
	    update_z_beta1_joint4(zg, beta, nalloc, n_acc, n_try, wtc, tau_eps, eta_up, eta_down, 
				  lambda_up, lambda_down, tau, gamma1, xx, indtau, ybar, ydata, like_choice, 
				  ngenes, nconds, nreps, ncomps, neffects, jstar, rand);
	    if(move_choice_eta==1) update_eta(eta_up, eta_down, lambda_up, lambda_down, aa_eta, bb_eta, 
					      zg, beta, nalloc, ngenes, nconds, jstar, rand);
	    if(move_choice_lam==1) update_lambda(lambda_up, lambda_down, eta_up, eta_down, 
						 lam1, lam2, nlam, zg, beta, nalloc, ngenes, jstar, rand);
	  }	  
	  update_wtc(wtc, nalloc, nu, ncomps, rand);
	}
	
	// these not part of MCMC chain so only do at thinned iterations
	if(iter%nthin==0){
	  if(like_choice==1) predict(ybar_pred1, ybar_pred2, ybar_pred3, ybar_pred4, ss_pred1, ss_pred2,
              pval_post_ss, pval_mix_ss, pval_post_ybar, pval_mix1_ybar, pval_mix2_ybar, 
              pval_mix3_ybar, pval_partial_ss, pval_partial_ybar, norm_ss, norm_ybar,
	      ybar, ss, tau, gamma1, aa, bb, zg, beta, xx, indtau, wtc, tau_eps, eta_up, eta_down, 
              lambda_up, lambda_down, 
              jstar, move_choice_bz, move_choice_tau, ngenes, nconds, nreps, neffects, rand);
	  else if(like_choice==2) predict_t(ybar_pred1, ybar_pred2, ybar_pred3, ybar_pred4, ss_pred1, 
	      ss_pred2, pval_post_ss, pval_mix_ss, pval_post_ybar, pval_mix1_ybar, pval_mix2_ybar, 
              pval_mix3_ybar, 
              ybar, ss, tau, gamma1, df, aa, bb, zg, beta, xx, indtau, wtc, tau_eps, eta_up, eta_down, 
              lambda_up, lambda_down, 
              jstar, move_choice_bz, move_choice_tau, ngenes, nconds, nreps, neffects, rand);
	  deviance_calc(deviance1, deviance2, beta, tau, gamma1, xx, indtau, ybar, ss, ydata, 
			df, like_choice, ngenes, nconds, nreps, neffects);
	}

	///////////// calc. means, probs and output traces	

	if(iter%nthin==0){
		//std::cerr << "out " << iter << "out ";
		for(int i=0; i<ngenes; ++i) {
		  mean_zg[i] += zg[i]-1;
		  if(zg[i]==0) prob0[i] += 1;
		  else if (zg[i]==2) prob2[i] += 1;
		  else prob1[i] += 1;
		  for(int j=0; j<neffects; ++j) {
			mean_beta[i][j] += beta[i][j];
		  }
		  for(int j=0; j<ntau; ++j) {
			mean_sig2[i][j] += 1.0/tau[i][j];
			mean_tau[i][j] += tau[i][j];
			meanlogscale_tau[i][j] += log(tau[i][j]);
		  }
		  if(like_choice==2){
		    for(int j=0; j<nrepstot; ++j) {
		      mean_gamma1[i][j] += gamma1[i][j];
		      meanlogscale_gamma1[i][j] += log(gamma1[i][j]);
		    }
		  }
		}
		for(int j=0; j<ntau; ++j) {
			mean_bb[j] += bb[j];
			if(move_choice_aa==1){
				mean_aa[j] += aa[j];
			}
		}
		for(int j=0; j<ncomps; ++j) {
			mean_wtc[j] += wtc[j];
		}
		mean_eta_up += eta_up;
		mean_eta_down += eta_down;
		mean_lambda_up += lambda_up;
		mean_lambda_down += lambda_down;
		mean_dev1 += deviance1;
		mean_dev2 += deviance2;

		/// trace files for model parameters
		if(trace_out==1){
		  for(int i=0; i<ngenes; ++i) {
		    zg_file << zg[i] << " ";
		    for(int j=0; j<neffects; ++j) {
		      beta_file << beta[i][j] << " ";
		    }
		    for(int j=0; j<ntau; ++j) {
		      sig2_file << 1.0/tau[i][j] << " ";
		    }
		  }
		  for(int j=0; j<ntau; ++j) {
		    bb_file << bb[j] << " ";
		    if(move_choice_aa==1){
		      aa_file << aa[j] << " ";
		    }
		  }
		  for(int j=0; j<ncomps; ++j) {
		    wtc_file << wtc[j] << " ";
		  }
		  eta_file << eta_up << " " << eta_down << " ";
		  if(move_choice_lam==1) lambda_file << lambda_up << " " << lambda_down << " ";
		  beta_file << std::endl;
		  sig2_file << std::endl;
		  bb_file << std::endl;
		}
		/// trace files for predicted data
		if(trace_pred==1){
		   for(int i=0; i<ngenes; ++i) {
		    for(int j=0; j<nconds; ++j) {
			ybar_pred1_file << ybar_pred1[i][j] << " ";
			ybar_pred2_file << ybar_pred2[i][j] << " ";
			ybar_pred3_file << ybar_pred3[i][j] << " ";
			ybar_pred4_file << ybar_pred4[i][j] << " ";
			ss_pred1_file << ss_pred1[i][j] << " ";
			ss_pred2_file << ss_pred2[i][j] << " ";
		    }}
		}
	}

} // end mcmc updates

//std::cerr << std::endl << "Number of samples recorded " << niter/nthin << std::endl;


//if(move_choice_bz!=1){
//	std::cout << "acceptance rate bz = " << n_acc << " " << n_try << " " << double(n_acc)/double(n_try) << std::endl;
//	summary_file << "acceptance rate bz = " << n_acc << " " << n_try << " " << double(n_acc)/double(n_try) << std::endl;
//}
//if(move_choice_aa==1){
//	std::cout << "acceptance rate a = " << n_aa_acc << " " << n_aa_try << " " << 
//                    double(n_aa_acc)/double(n_aa_try) << std::endl;
//	summary_file << "acceptance rate a = " << n_aa_acc << " " << n_aa_try << " " << 
//                       double(n_aa_acc)/double(n_aa_try) << std::endl;
//}
//if(move_choice_tau==2){
//	std::cout << "acceptance rate tau logNorm = " << n_tau_lognorm_acc << " " << 
//	  n_tau_lognorm_try << " " << 
//	  double(n_tau_lognorm_acc)/double(n_tau_lognorm_try) << std::endl;
//	summary_file << "acceptance rate tau logNorm = " << n_tau_lognorm_acc << " " << 
//	  n_tau_lognorm_try << " " << 
//	  double(n_tau_lognorm_acc)/double(n_tau_lognorm_try) << std::endl;
//}

// std::cout << run_dir << " Output Directory" << std::endl;
// PRINTF("\n Output Directory is %s \n",dirname);
// PRINTF("\n Output Directory \n");
PRINTF("\n");


////////////////////////////// divide means, probs by no. samples ////////////////////////////////

/// pvals are thinned now
for(int i=0; i<ngenes; ++i) {
for(int j=0; j<nconds; ++j) {
	pval_post_ss[i][j] = pval_post_ss[i][j] / float(niter/nthin);
	pval_mix_ss[i][j] = pval_mix_ss[i][j] / float(niter/nthin);
}}
////// for predictive p-vals conditioned on zg, must divide by prob(zg=0) etc.
for(int i=0; i<ngenes; ++i) {
	pval_post_ybar[i][0] = pval_post_ybar[i][0] / (prob0[i]);
	pval_mix1_ybar[i][0] = pval_mix1_ybar[i][0] / (prob0[i]);
	pval_mix2_ybar[i][0] = pval_mix2_ybar[i][0] / (prob0[i]);
	pval_mix3_ybar[i][0] = pval_mix3_ybar[i][0] / (prob0[i]);
	pval_post_ybar[i][1] = pval_post_ybar[i][1] / (prob1[i]);
	pval_mix1_ybar[i][1] = pval_mix1_ybar[i][1] / (prob1[i]);
	pval_mix2_ybar[i][1] = pval_mix2_ybar[i][1] / (prob1[i]);
	pval_mix3_ybar[i][1] = pval_mix3_ybar[i][1] / (prob1[i]);
	pval_post_ybar[i][2] = pval_post_ybar[i][2] / (prob2[i]);
	pval_mix1_ybar[i][2] = pval_mix1_ybar[i][2] / (prob2[i]);
	pval_mix2_ybar[i][2] = pval_mix2_ybar[i][2] / (prob2[i]);
	pval_mix3_ybar[i][2] = pval_mix3_ybar[i][2] / (prob2[i]);
}
//// partial p-vals: divide by norm
if(like_choice==1){
  for(int i=0; i<ngenes; ++i) {
    for(int j=0; j<nconds; ++j) {
      pval_partial_ss[i][j] = pval_partial_ss[i][j] / norm_ss[i][j];
    }
    pval_partial_ybar[i] = pval_partial_ybar[i] / norm_ybar[i];
  }
}

for(int i=0; i<ngenes; ++i) {
  mean_zg[i] = mean_zg[i]/float(niter/nthin);
  prob0[i] = prob0[i]/float(niter/nthin);
  prob1[i] = prob1[i]/float(niter/nthin);
  prob2[i] = prob2[i]/float(niter/nthin);
  for(int j=0; j<neffects; ++j) {
	mean_beta[i][j] = mean_beta[i][j]/float(niter/nthin);
}}

for(int i=0; i<ngenes; ++i) {
  for(int j=0; j<ntau; ++j) {
    mean_sig2[i][j] = mean_sig2[i][j]/float(niter/nthin);
    mean_tau[i][j] = mean_tau[i][j]/float(niter/nthin);
    meanlogscale_tau[i][j] = exp(meanlogscale_tau[i][j]/float(niter/nthin));
  }
  if(like_choice==2){
    for(int j=0; j<nrepstot; ++j) {
      mean_gamma1[i][j] = mean_gamma1[i][j]/float(niter/nthin);
      meanlogscale_gamma1[i][j] = exp(meanlogscale_gamma1[i][j]/float(niter/nthin));
    }
  }
}

for(int j=0; j<ntau; ++j) mean_bb[j] = mean_bb[j]/float(niter/nthin);
for(int j=0; j<ntau; ++j) mean_aa[j] = mean_aa[j]/float(niter/nthin);
for(int j=0; j<ncomps; ++j) mean_wtc[j] = mean_wtc[j]/float(niter/nthin);

mean_eta_up = mean_eta_up/float(niter/nthin);
mean_eta_down = mean_eta_down/float(niter/nthin);
mean_lambda_up = mean_lambda_up/float(niter/nthin);
mean_lambda_down = mean_lambda_down/float(niter/nthin);

mean_dev1 = mean_dev1/float(niter/nthin);
mean_dev2 = mean_dev2/float(niter/nthin);
double dthetahat1 = 0, dthetahat2 = 0;
double dum1,dum2;
 int ind;
for(int g=0; g<ngenes; ++g) {
for(int c=0; c<nconds; ++c) {
  dum2 = 0;
  for(int j=0; j<neffects; ++j) {
    dum2 += mean_beta[g][j]*xx[j][c];
  }
  if(like_choice==1){
    dum1 = (nreps[c]-1)*ss[g][c] + nreps[c]*pow(ybar[g][c] - dum2, 2);
    dthetahat1 += meanlogscale_tau[g][indtau[c]]*dum1 - nreps[c]*log(meanlogscale_tau[g][indtau[c]]);
  }
  else if(like_choice==2){
    for(int r=0; r<nreps[c]; ++r) {
      if(c==0) ind = r;
      else ind = c*nreps[c-1] + r;
      dum1 = pow(ydata[g][ind] - dum2, 2);
      dthetahat1 += meanlogscale_tau[g][indtau[c]]*meanlogscale_gamma1[g][ind]*dum1 - 
	log(meanlogscale_gamma1[g][ind]*meanlogscale_tau[g][indtau[c]]);
      dthetahat2 += (df[c]+1)*log(1 + dum1*meanlogscale_tau[g][indtau[c]]/df[c]) - 
	log(df[c]*meanlogscale_tau[g][indtau[c]]);
    }
  }
}}

///////////////////////////////////////////// output ///////////////////////////////////

for(int i=0; i<ngenes; ++i) {
  mzg_file << mean_zg[i] << " ";
  probclass_file << prob0[i] << " " << prob1[i] << " " << prob2[i] << " ";
  for(int j=0; j<neffects; ++j) {
	mbeta_file << mean_beta[i][j] << " ";
}}

for(int i=0; i<ngenes; ++i) {
  for(int j=0; j<ntau; ++j) {
	msig2_file << mean_sig2[i][j] << " ";
	mtau_file << mean_tau[i][j] << " ";
  }
  for(int j=0; j<nconds; ++j) {
	post_ss_file << pval_post_ss[i][j] << " ";
	mix_ss_file << pval_mix_ss[i][j] << " ";
	if(like_choice==1) part_ss_file << pval_partial_ss[i][j] << " ";
  }
}

for(int i=0; i<ngenes; ++i) {
for(int j=0; j<ncomps; ++j) {
	post_ybar_file << pval_post_ybar[i][j] << " ";
	mix1_ybar_file << pval_mix1_ybar[i][j] << " ";
	mix2_ybar_file << pval_mix2_ybar[i][j] << " ";
	mix3_ybar_file << pval_mix3_ybar[i][j] << " ";
}
if(like_choice==1) part_ybar_file << pval_partial_ybar[i] << " ";
}

for(int j=0; j<ntau; ++j) mbb_file << mean_bb[j] << " ";
if(move_choice_aa==1) for(int j=0; j<ntau; ++j) maa_file << mean_aa[j] << " ";
for(int j=0; j<ncomps; ++j) mwtc_file << mean_wtc[j] << " ";

m_eta_file << mean_eta_up << " " << mean_eta_down << std::endl;
m_lambda_file << mean_lambda_up << " " << mean_lambda_down << std::endl;


summary_file << "post. mean deviance (Dbar) = " << mean_dev1 << std::endl;
summary_file << "D(thetahat) = " << dthetahat1 << std::endl;
summary_file << "p_D = " << mean_dev1 - dthetahat1 << std::endl;
summary_file << "DIC = " << 2*mean_dev1 - dthetahat1 << std::endl;
summary_file << "post. mean deviance (Dbar) (2) = " << mean_dev2 << std::endl;
summary_file << "D(thetahat) (2) = " << dthetahat2 << std::endl;
summary_file << "p_D (2) = " << mean_dev2 - dthetahat2 << std::endl;
summary_file << "DIC (2) = " << 2*mean_dev2 - dthetahat2 << std::endl;


summary_file.close();

beta_file.close();
sig2_file.close();
bb_file.close();
aa_file.close();
wtc_file.close();
zg_file.close();
eta_file.close();
lambda_file.close();

ybar_pred1_file.close();
ybar_pred2_file.close();
ybar_pred3_file.close();
ybar_pred4_file.close();
ss_pred1_file.close();
ss_pred2_file.close();

mbeta_file.close();
msig2_file.close();
mbb_file.close();
maa_file.close();
mtau_file.close();
mwtc_file.close();
mzg_file.close();
m_eta_file.close();
m_lambda_file.close();

probclass_file.close();

post_ss_file.close();
mix_ss_file.close();
part_ss_file.close();
post_ybar_file.close();
mix1_ybar_file.close();
mix2_ybar_file.close();
mix3_ybar_file.close();
part_ybar_file.close();

}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

void freeBGmixMemory(int &ngenes, int &neffects){

  for(int i=0; i<ngenes; ++i) delete [] ybar[i];
  delete [] ybar;
  for(int i=0; i<ngenes; ++i) delete [] ss[i];
  delete [] ss;
  for(int i=0; i<ngenes; ++i) delete [] ydata[i];
  delete [] ydata;

  for(int i=0; i<ngenes; ++i) delete [] beta[i];
  delete [] beta;
  for(int i=0; i<ngenes; ++i) delete [] tau[i];
  delete [] tau;
  for(int i=0; i<ngenes; ++i) delete [] gamma1[i];
  delete [] gamma1;

  delete [] zg;

  for(int i=0; i<ngenes; ++i) delete [] mean_beta[i];
  delete [] mean_beta;
  for(int i=0; i<ngenes; ++i) delete [] mean_tau[i];
  delete [] mean_tau;
  for(int i=0; i<ngenes; ++i) delete [] meanlogscale_tau[i];
  delete [] meanlogscale_tau;
  for(int i=0; i<ngenes; ++i) delete [] mean_sig2[i];
  delete [] mean_sig2;
  for(int i=0; i<ngenes; ++i) delete [] mean_gamma1[i];
  delete [] mean_gamma1;
  for(int i=0; i<ngenes; ++i) delete [] meanlogscale_gamma1[i];
  delete [] meanlogscale_gamma1;

  delete [] mean_zg;
  delete [] prob0;
  delete [] prob1;
  delete [] prob2;

  for(int i=0; i<ngenes; ++i) delete [] ybar_pred1[i];
  delete [] ybar_pred1;
  for(int i=0; i<ngenes; ++i) delete [] ybar_pred2[i];
  delete [] ybar_pred2;
  for(int i=0; i<ngenes; ++i) delete [] ybar_pred3[i];
  delete [] ybar_pred3;
  for(int i=0; i<ngenes; ++i) delete [] ybar_pred4[i];
  delete [] ybar_pred4;
  for(int i=0; i<ngenes; ++i) delete [] ss_pred1[i];
  delete [] ss_pred1;
  for(int i=0; i<ngenes; ++i) delete [] ss_pred2[i];
  delete [] ss_pred2;

  for(int i=0; i<ngenes; ++i) delete [] pval_post_ss[i];
  delete [] pval_post_ss;
  for(int i=0; i<ngenes; ++i) delete [] pval_mix_ss[i];
  delete [] pval_mix_ss;
  for(int i=0; i<ngenes; ++i) delete [] pval_partial_ss[i];
  delete [] pval_partial_ss;

  for(int i=0; i<ngenes; ++i) delete [] pval_post_ybar[i];
  delete [] pval_post_ybar;
  for(int i=0; i<ngenes; ++i) delete [] pval_mix1_ybar[i];
  delete [] pval_mix1_ybar;
  for(int i=0; i<ngenes; ++i) delete [] pval_mix2_ybar[i];
  delete [] pval_mix2_ybar;
  for(int i=0; i<ngenes; ++i) delete [] pval_mix3_ybar[i];
  delete [] pval_mix3_ybar;
  for(int i=0; i<ngenes; ++i) delete [] norm_ss[i];
  delete [] norm_ss;

  delete [] pval_partial_ybar;
  delete [] norm_ybar;

  delete [] aa;
  delete [] bb;
  delete [] nu;
  delete [] wtc;
  delete [] df;
  delete [] nalloc;

  for(int i=0; i<neffects; ++i) delete [] xx[i];
  delete [] xx;

  delete [] mean_aa;
  delete [] mean_bb;
  delete [] mean_wtc;

}






