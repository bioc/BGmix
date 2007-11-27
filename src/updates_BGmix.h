#include "rand.hh"

void update_beta0(double** beta, double** tau, double** gamma,
		  double** xx, int* indtau, double** ybar, double** ydata, 
		  int &like_choice, int &ngenes, 
		  int &nconds, int* nreps, int &neffects, 
		  int &jstar, Random &rand, std::ofstream &summary_file);

void update_z_beta1_joint1(int* zg, double* wtc, int* nalloc, double &eta_up, double &eta_down, 
			   double** beta, double** tau, double** gamma, double** xx, 
			   int* indtau, double** ybar, double** ydata, int &like_choice, int &ngenes, 
			   int &nconds, int* nreps, 
			   int &ncomps, int &neffects, int &jstar, Random &rand);

void update_z_beta1_joint2(int* zg, double** beta, int* nalloc, int &n_acc, int &n_try, double* wtc, 
			   double &eta_up, double &eta_down, double** tau, double** gamma,
			   double** xx, int* indtau, double** ybar, double** ydata, int &like_choice, 
			   int &ngenes, int &nconds, int* nreps, int &ncomps, int &neffects, 
			   int &jstar, Random &rand);

void update_z_beta1_joint3(int* zg, double** beta, int* nalloc, int &n_acc, int &n_try, double* wtc, 
			   double &eta_up, double &eta_down, double &lambda_up, double &lambda_down, 
			   double** tau, double** gamma, double** xx, int* indtau, 
			   double** ybar, double** ydata, int &like_choice, int &ngenes, 
			   int &nconds, int* nreps, int &ncomps, int &neffects, int &jstar, Random &rand);

void update_z_beta1_joint4(int* zg, double** beta, int* nalloc, int &n_acc, int &n_try, 
			   double* wtc, double &tau_eps, double &eta_up, double &eta_down, 
			   double &lambda_up, double &lambda_down, double** tau, double** gamma, 
			   double** xx, int* indtau, double** ybar, double** ydata, int &like_choice, 
			   int &ngenes, int &nconds, int* nreps, 
			   int &ncomps, int &neffects, int &jstar, Random &rand);

void update_tau(double** beta, double** tau, double** gamma,
		double** xx, int* indtau, double** ybar, double** ss, double** ydata, 
		double* aa, double* bb, int &like_choice, int &ngenes, int &nconds, 
		int &ntau, int* nreps, int &neffects, Random &rand);

void update_tau_logNorm(double** beta, double** tau, double** gamma,
			double** xx, int* indtau, double** ybar, double** ss, double** ydata, 
			double* aa, double* bb, 
			int &n_tau_lognorm_acc, int &n_tau_lognorm_try,
			int &like_choice, int &ngenes, int &nconds, int &ntau, int* nreps, 
			int &neffects, Random &rand);
  
void update_tau_cut(double** tau, double** ss, int* indtau, double* aa, double* bb, 
		    int &ngenes, int &nconds, int &ntau, int* nreps, Random &rand);

void update_gamma(double** beta, double** tau, double** gamma,
		  double** xx, int* indtau, double** ydata, double* df, 
		  int &ngenes, int &nconds, int &ntau, int* nreps, 
		  int &neffects, Random &rand);


/////////////////

void update_bb(double** tau, double* aa, double* bb, double &gg, double &hh, 
	       int &ngenes, int &nconds, int &ntau, Random &rand);

void update_bb_logNorm(double** tau, double* aa, double* bb, double &gg, double &hh, 
		       int &ngenes, int &nconds, int &ntau, Random &rand);
  
void update_aa(double &sig_aa, double** tau, double* aa, double* bb, double &gg, double &hh, 
	       int &n_aa_acc, int &n_aa_try, int &ngenes, int &nconds, int &ntau, Random &rand);

void update_aa_logNorm(double &sig_aa, double** tau, double* aa, double* bb, double &tau_var, 
		       int &ngenes, int &nconds, int &ntau, Random &rand);


void update_eta_unif(double &eta_up, double &eta_down, double &aa_eta, double &bb_eta, 
		     int* nalloc, Random &rand);

void update_eta(double &eta_up, double &eta_down, double &lambda_up, double &lambda_down, 
		double &aa_eta, double &bb_eta, int* zg, double** beta, int* nalloc, 
		int &ngenes, int &nconds, int &jstar, Random &rand);

void update_lambda(double &lambda_up, double &lambda_down, double &eta_up, double &eta_down, 
		   double &lam1, double &lam2, int &nlam, int* zg, double** beta, int* nalloc, 
		   int &ngenes, int &jstar, Random &rand);

void update_wtc(double* wtc, int* nalloc, double* nu, int &ncomps, Random &rand);


/////////////////

void predict(double** ybar_pred1, double** ybar_pred2, double** ybar_pred3, double** ybar_pred4, 
double** ss_pred1, double** ss_pred2, double** pval_post_ss, double** pval_mix_ss, 
double** pval_post_ybar, double** pval_mix1_ybar, double** pval_mix2_ybar, double** pval_mix3_ybar, 
//double* pval_post_ybar, double* pval_mix1_ybar, double* pval_mix2_ybar, double* pval_mix3_ybar, 
double** pval_partial_ss, double* pval_partial_ybar, double** norm_ss, double* norm_ybar,
double** ybar, double** ss, 
double** tau, double** gamma, double* aa, double* bb, 
int* zg, double** beta, double** xx, int* indtau, double* wtc, double &tau_eps,
double &eta_up, double &eta_down, double &lambda_up, double &lambda_down, int &jstar, 
int &move_choice_bz, int &move_choice_tau, int &ngenes, 
int &nconds, int* nreps, int &neffects, Random &rand);

void predict_t(double** ybar_pred1, double** ybar_pred2, double** ybar_pred3, double** ybar_pred4, 
double** ss_pred1, double** ss_pred2, double** pval_post_ss, double** pval_mix_ss, 
double** pval_post_ybar, double** pval_mix1_ybar, double** pval_mix2_ybar, double** pval_mix3_ybar, 
double** ybar, double** ss, double** tau, double** gamma, double* df, double* aa, double* bb, 
int* zg, double** beta, double** xx, int* indtau, double* wtc, double &tau_eps,
double &eta_up, double &eta_down, double &lambda_up, double &lambda_down, int &jstar, 
int &move_choice_bz, int &move_choice_tau, int &ngenes, 
	       int &nconds, int* nreps, int &neffects, Random &rand);

void deviance_calc(double &deviance1, double &deviance2, double** beta, double** tau, 
		   double** gamma, double** xx, int* indtau, double** ybar, 
		   double** ss, double** ydata, double* df, int &like_choice, int &ngenes, 
		   int &nconds, int* nreps, int &neffects);










