using namespace std;

#ifndef _BGMIX_H
#define _BGMIX_H
extern "C"
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
	   char** dirname, char **basepath);

extern "C"
void freeBGmixMemory(int &ngenes, int &neffects);

#endif // _BGMIX_H




