#include <cstdlib>
#include <iostream>
#include "updates_BGmix.h"
#include "rand.hh"
#include <Rmath.h>
#include <valarray>

/////////////////////////////////////////////////////////////////////////
//// Gibbs for beta0 (flat prior)

//// arrays start at 0

void update_beta0(double** beta, double** tau, double** gamma,
		  double** xx, int* indtau, double** ybar, double** ydata, 
		  int &like_choice, int &ngenes, 
		  int &nconds, int* nreps, int &neffects, 
		  int &jstar, Random &rand, std::ofstream &summary_file){

double dum2,dum3,dum4;
int ind;

//std::cerr << "beta0" << " ";

for(int j=0; j<neffects; ++j) {
if(j!=jstar){
for(int g=0; g<ngenes; ++g) {
  dum3 = 0;
  dum4 = 0;
  for(int c=0; c<nconds; ++c) {
    dum2 = 0;
    for(int jp=0; jp<neffects; ++jp) {
      if(jp!=j) dum2 += beta[g][jp]*xx[jp][c];
    }
    if(like_choice==1){
      dum3 += tau[g][indtau[c]]*nreps[c]*xx[j][c]*(ybar[g][c]-dum2);
      dum4 += tau[g][indtau[c]]*nreps[c]*pow(xx[j][c],2);
    }
    else if(like_choice==2){
      for(int r=0; r<nreps[c]; ++r) {
	if(c==0) ind = r;
	else ind = c*nreps[c-1] + r;
	dum3 += tau[g][indtau[c]]*xx[j][c]*gamma[g][ind]*(ydata[g][ind]-dum2);
	dum4 += tau[g][indtau[c]]*pow(xx[j][c],2)*gamma[g][ind];
      }
    }
  }
  beta[g][j] = rand.Normal( dum3/dum4 , sqrt(1.0/(dum4))  );
}}}

// std::cerr << "  beta0 " << beta[0][0] << std::endl;

}

/////////////////////////////////////////////////////////////////////////
//// Joint Gibbs for alloc and beta1 (beta1 ~ Unif) 
//// with beta1 integrated out for alloc update

void update_z_beta1_joint1(int* zg, double* wtc, int* nalloc, double &eta_up, double &eta_down, 
			   double** beta, double** tau, double** gamma, double** xx, 
			   int* indtau, double** ybar, double** ydata, int &like_choice, int &ngenes, 
			   int &nconds, int* nreps, 
			   int &ncomps, int &neffects, int &jstar, Random &rand){

double dum2,dum3,dum4,palloc1,palloc2,palloc0,psum,uu,u_up,u_down;
int ind;

for(int k=0; k<ncomps; ++k) nalloc[k] = 0;

for(int g=0; g<ngenes; ++g) {

  dum3 = 0;
  dum4 = 0;
  for(int c=0; c<nconds; ++c) {
    dum2 = 0;
    for(int jp=0; jp<neffects; ++jp) {
      if(jp!=jstar) dum2 += beta[g][jp]*xx[jp][c];
    }
    if(like_choice==1){
      dum3 += tau[g][indtau[c]]*nreps[c]*xx[jstar][c]*(ybar[g][c]-dum2);
      dum4 += tau[g][indtau[c]]*nreps[c]*pow(xx[jstar][c],2);
    }
    else if(like_choice==2){
      for(int r=0; r<nreps[c]; ++r) {
	if(c==0) ind = r;
	else ind = c*nreps[c-1] + r;
	dum3 += tau[g][indtau[c]]*xx[jstar][c]*gamma[g][ind]*(ydata[g][ind]-dum2);
	dum4 += tau[g][indtau[c]]*pow(xx[jstar][c],2)*gamma[g][ind];
      }
    }
  }

  /////////// replaced gsl with R fns
    // palloc0 = gsl_sf_erf_Q( sqrt(1.0/dum4)*dum3 ) - gsl_sf_erf_Q( sqrt(1.0/dum4)*dum3 + 
    //							      eta_down*sqrt(dum4) );
    // palloc2 = gsl_sf_erf_Q( sqrt(1.0/dum4)*dum3 ) - gsl_sf_erf_Q( sqrt(1.0/dum4)*dum3 - 
    //								      eta_up*sqrt(dum4) );
    palloc0 = pnorm( sqrt(1.0/dum4)*dum3 ,0,1,0,0) - 
      pnorm( sqrt(1.0/dum4)*dum3 + eta_down*sqrt(dum4) ,0,1,0,0);
    palloc2 = pnorm( sqrt(1.0/dum4)*dum3 ,0,1,0,0) - 
      pnorm( sqrt(1.0/dum4)*dum3 - eta_up*sqrt(dum4) ,0,1,0,0);

	palloc0 = palloc0 * wtc[0] * sqrt(2.0*3.14159/(dum4)) / eta_down;
	palloc2 = -palloc2 * wtc[2] * sqrt(2.0*3.14159/(dum4)) / eta_up;
	palloc1 = wtc[1] * exp(-pow(dum3,2)/(dum4*2));
	psum = palloc0 + palloc1 + palloc2;

	uu = rand.Uniform()*psum;
	if(uu<palloc0) {
		zg[g] = 0;
		nalloc[0]++;
	}
	else if(uu<palloc0+palloc2) {
		zg[g] = 2;
		nalloc[2]++;
	}
	else {
		zg[g] = 1;
		nalloc[1]++;
	}

  /////////// replaced gsl with R fns
	if(zg[g]==1) beta[g][jstar] = 0;
	else if(zg[g]==0){
	  // u_down = gsl_cdf_gaussian_P( -eta_down-dum3/dum4 , sqrt(1.0/(dum4)) );
	  // u_up = gsl_cdf_gaussian_P( -dum3/dum4 , sqrt(1.0/(dum4)) );
	  u_down = pnorm( -eta_down-dum3/dum4 , 0, sqrt(1.0/(dum4)) ,1,0);
	  u_up = pnorm( -dum3/dum4 , 0, sqrt(1.0/(dum4)) ,1,0);
	  uu = rand.Uniform(u_down,u_up);
	  // beta[g][jstar] = gsl_cdf_gaussian_Pinv(uu, sqrt(1.0/(dum4)) ) + dum3/dum4; 
	  beta[g][jstar] = qnorm(uu, 0, sqrt(1.0/(dum4)) ,1,0) + dum3/dum4; 
	}
	else{
	  // u_down = gsl_cdf_gaussian_P( -dum3/dum4 , sqrt(1.0/(dum4)) );
	  // u_up = gsl_cdf_gaussian_P( eta_up-dum3/dum4 , sqrt(1.0/(dum4)) );
	  u_down = pnorm( -dum3/dum4 , 0, sqrt(1.0/(dum4)) ,1,0);
	  u_up = pnorm( eta_up-dum3/dum4 , 0, sqrt(1.0/(dum4)) ,1,0);
	  uu = rand.Uniform(u_down,u_up);
	  // beta[g][jstar] = gsl_cdf_gaussian_Pinv(uu, sqrt(1.0/(dum4)) ) + dum3/dum4; 
	  beta[g][jstar] = qnorm(uu, 0, sqrt(1.0/(dum4)) ,1,0) + dum3/dum4; 
	}
}

// std::cerr << "  beta1 " << beta[0][1] << std::endl;

}

/////////////////////////////////////////////////////////////////////////
//// Joint MH for alloc and beta1 (beta1 ~ Unif) 

void update_z_beta1_joint2(int* zg, double** beta, int* nalloc, int &n_acc, int &n_try, double* wtc, 
			   double &eta_up, double &eta_down, double** tau, double** gamma,
			   double** xx, int* indtau, double** ybar, double** ydata, int &like_choice, 
			   int &ngenes, int &nconds, int* nreps, int &ncomps, int &neffects, 
			   int &jstar, Random &rand){

//std::cerr << "joint2" << " ";

double dum2,dum3,dum4,phi,beta1try,uu,acc_prob;
int zg_try,indic;
int ind;

for(int g=0; g<ngenes; ++g) {

  dum3 = 0;
  dum4 = 0;
  for(int c=0; c<nconds; ++c) {
    dum2 = 0;
    for(int jp=0; jp<neffects; ++jp) {
      if(jp!=jstar) dum2 += beta[g][jp]*xx[jp][c];
    }
    if(like_choice==1){
      dum3 += tau[g][indtau[c]]*nreps[c]*xx[jstar][c]*(ybar[g][c]-dum2);
      dum4 += tau[g][indtau[c]]*nreps[c]*pow(xx[jstar][c],2);
    }
    else if(like_choice==2){
      for(int r=0; r<nreps[c]; ++r) {
	if(c==0) ind = r;
	else ind = c*nreps[c-1] + r;
	dum3 += tau[g][indtau[c]]*xx[jstar][c]*gamma[g][ind]*(ydata[g][ind]-dum2);
	dum4 += tau[g][indtau[c]]*pow(xx[jstar][c],2)*gamma[g][ind];
      }
    }
  }

	phi = exp(-0.5*dum4*pow(dum3/dum4,2)) * sqrt(dum4/6.28318);
	indic = 1;

	uu = rand.Uniform();
	if(uu<wtc[0]) {
		zg_try = 0;
		beta1try = rand.Normal( dum3/dum4 , sqrt(1.0/(dum4))  );
		if(beta1try>0 | beta1try<-eta_down) indic = 0;
		if(zg[g]==0) acc_prob = indic;
		else if(zg[g]==2) acc_prob = indic*eta_up/eta_down; 
		else acc_prob = indic/(phi*eta_down);
	}
	else if(uu<wtc[0]+wtc[2]) {
		zg_try = 2;
		beta1try = rand.Normal( dum3/dum4 , sqrt(1.0/(dum4))  );
		if(beta1try<0 | beta1try>eta_up) indic = 0;
		if(zg[g]==0) acc_prob = indic*eta_down/eta_up; 
		else if(zg[g]==2) acc_prob = indic;
		else acc_prob = indic/(phi*eta_up);
	}
	else {
		zg_try = 1;
		beta1try = 0;
		if(zg[g]==0) acc_prob = phi*eta_down;
		else if(zg[g]==2) acc_prob = phi*eta_up;
		else acc_prob = 1;
	}

//std::cerr << zg[g] << " " << zg_try << " " << acc_prob << "    ";

	uu = rand.Uniform();
	n_try++;
	if(uu<acc_prob){
		nalloc[zg_try]++;
		nalloc[zg[g]]--;
		zg[g] = zg_try;
		beta[g][jstar] = beta1try;
		n_acc++; 
	}

}

//std::cerr << "sum alloc" << nalloc[0]+nalloc[1]+nalloc[2] << " ";

}

/////////////////////////////////////////////////////////////////////////
//// Joint MH for alloc and beta1 (beta1 ~ Gamma) 

void update_z_beta1_joint3(int* zg, double** beta, int* nalloc, int &n_acc, int &n_try, double* wtc, 
			   double &eta_up, double &eta_down, double &lambda_up, double &lambda_down, 
			   double** tau, double** gamma, double** xx, int* indtau, 
			   double** ybar, double** ydata, int &like_choice, int &ngenes, 
			   int &nconds, int* nreps, int &ncomps, int &neffects, int &jstar, Random &rand){

//std::cerr << "joint3" << " ";

double dum2,dum3,dum4,phi,phigam_up,phigam_down,ratiogam,beta1try,uu,acc_prob;
int zg_try;
int ind;

for(int g=0; g<ngenes; ++g) {

  dum3 = 0;
  dum4 = 0;
  for(int c=0; c<nconds; ++c) {
    dum2 = 0;
    for(int jp=0; jp<neffects; ++jp) {
      if(jp!=jstar) dum2 += beta[g][jp]*xx[jp][c];
    }
    if(like_choice==1){
      dum3 += tau[g][indtau[c]]*nreps[c]*xx[jstar][c]*(ybar[g][c]-dum2);
      dum4 += tau[g][indtau[c]]*nreps[c]*pow(xx[jstar][c],2);
    }
    else if(like_choice==2){
      for(int r=0; r<nreps[c]; ++r) {
	if(c==0) ind = r;
	else ind = c*nreps[c-1] + r;
	dum3 += tau[g][indtau[c]]*xx[jstar][c]*gamma[g][ind]*(ydata[g][ind]-dum2);
	dum4 += tau[g][indtau[c]]*pow(xx[jstar][c],2)*gamma[g][ind];
      }
    }
  }

	phi = exp(-0.5*dum4*pow(dum3/dum4,2)) * sqrt(dum4/6.28318);
  /////////// replaced gsl with R fns
	  // phigam_up = phi*gsl_sf_gamma(lambda_up);
	  // phigam_down = phi*gsl_sf_gamma(lambda_down);
	phigam_up = phi*gammafn(lambda_up);
	phigam_down = phi*gammafn(lambda_down);
	ratiogam = phigam_up/phigam_down;

	uu = rand.Uniform();
	if(uu<wtc[0]) {
	  zg_try = 0;
	  beta1try = rand.Normal( dum3/dum4 , sqrt(1.0/(dum4))  );
	  if(beta1try>0) acc_prob = 0;
	  else{
	    if(zg[g]==0) acc_prob = pow(beta1try/beta[g][jstar],lambda_down-1) * 
	      exp( eta_down*(beta1try-beta[g][jstar]) );
	    else if(zg[g]==2) acc_prob = (eta_down/eta_up) * 
	      exp(eta_down*beta1try+eta_up*beta[g][jstar]) *
	      ratiogam * pow(-beta1try*eta_down,lambda_down-1)/pow(beta[g][jstar]*eta_up,lambda_up-1);
	    else acc_prob = pow(-beta1try*eta_down,lambda_down-1) * eta_down *
	      exp(eta_down*beta1try)/phigam_down;
	  }
	}
	else if(uu<wtc[0]+wtc[2]) {
	  zg_try = 2;
	  beta1try = rand.Normal( dum3/dum4 , sqrt(1.0/(dum4))  );
	  if(beta1try<0) acc_prob = 0;
	  else{	    
	    if(zg[g]==0) acc_prob = (eta_up/eta_down) * exp(-eta_up*beta1try-eta_down*beta[g][jstar]) *
	      (1./ratiogam) * pow(beta1try*eta_up,lambda_up-1)/pow(-beta[g][jstar]*eta_down,lambda_down-1);
	    else if(zg[g]==2) acc_prob = pow(beta1try/beta[g][jstar],lambda_up-1) * 
	      exp(-eta_up*(beta1try-beta[g][jstar]));
	    else acc_prob = pow(beta1try*eta_up,lambda_up-1) * eta_up *
	      exp(-eta_up*beta1try)/phigam_up;
	  }
	}
	else {
	  zg_try = 1;
	  beta1try = 0;
	  if(zg[g]==0) acc_prob = exp(-eta_down*beta[g][jstar])*phigam_down / 
	    ( eta_down * pow(-beta[g][jstar]*eta_down,lambda_down-1) );
	  else if(zg[g]==2) acc_prob = exp(eta_up*beta[g][jstar])*phigam_up / 
	    ( eta_up * pow(beta[g][jstar]*eta_up,lambda_up-1) );
	  else acc_prob = 1;
	}

	//	if(acc_prob<0) std::cerr << "  z " << zg[g] << " " << zg_try << "   beta " << beta[g][jstar] << 
	//	  " " << beta1try << " acc prob " << acc_prob << std::endl;

	uu = rand.Uniform();
	n_try++;
	if(uu<acc_prob){
		nalloc[zg_try]++;
		nalloc[zg[g]]--;
		zg[g] = zg_try;
		beta[g][jstar] = beta1try;
		n_acc++; 
	}

}

//std::cerr << "sum alloc" << nalloc[0]+nalloc[1]+nalloc[2] << " ";

}

/////////////////////////////////////////////////////////////////////////
//// Joint MH for alloc and beta1 (beta1 ~ Gamma) 
//// This one has N(0,epsilon) for null

void update_z_beta1_joint4(int* zg, double** beta, int* nalloc, int &n_acc, int &n_try, 
			   double* wtc, double &tau_eps, double &eta_up, double &eta_down, 
			   double &lambda_up, double &lambda_down, double** tau, double** gamma, 
			   double** xx, int* indtau, double** ybar, double** ydata, int &like_choice, 
			   int &ngenes, int &nconds, int* nreps, 
			   int &ncomps, int &neffects, int &jstar, Random &rand){

//std::cerr << "joint4" << " ";

double dum2,dum3,dum4,gam_up,gam_down,ratiogam,beta1try,uu,acc_prob;
int zg_try;
int ind;

for(int g=0; g<ngenes; ++g) {

  dum3 = 0;
  dum4 = 0;
  for(int c=0; c<nconds; ++c) {
    dum2 = 0;
    for(int jp=0; jp<neffects; ++jp) {
      if(jp!=jstar) dum2 += beta[g][jp]*xx[jp][c];
    }
    if(like_choice==1){
      dum3 += tau[g][indtau[c]]*nreps[c]*xx[jstar][c]*(ybar[g][c]-dum2);
      dum4 += tau[g][indtau[c]]*nreps[c]*pow(xx[jstar][c],2);
    }
    else if(like_choice==2){
      for(int r=0; r<nreps[c]; ++r) {
	if(c==0) ind = r;
	else ind = c*nreps[c-1] + r;
	dum3 += tau[g][indtau[c]]*xx[jstar][c]*gamma[g][ind]*(ydata[g][ind]-dum2);
	dum4 += tau[g][indtau[c]]*pow(xx[jstar][c],2)*gamma[g][ind];
      }
    }
  }

       	beta1try = rand.Normal( dum3/dum4 , sqrt(1.0/(dum4))  );

  /////////// replaced gsl with R fns
	  // gam_up = gsl_sf_gamma(lambda_up);
	  // gam_down = gsl_sf_gamma(lambda_down);
	gam_up = gammafn(lambda_up);
	gam_down = gammafn(lambda_down);
	ratiogam = gam_up/gam_down;

	uu = rand.Uniform();
	if(uu<wtc[0]) {
	  zg_try = 0;
	  if(beta1try>0) acc_prob = 0;
	  else{
	    if(zg[g]==0) acc_prob = pow(beta1try/beta[g][jstar],lambda_down-1) * 
	      exp( eta_down*(beta1try-beta[g][jstar]) );
	    else if(zg[g]==2) acc_prob = (eta_down/eta_up) * 
	      exp(eta_down*beta1try+eta_up*beta[g][jstar]) * ratiogam * 
	      pow(-beta1try*eta_down,lambda_down-1)/pow(beta[g][jstar]*eta_up,lambda_up-1);
	    else acc_prob = pow(-beta1try*eta_down,lambda_down-1) * eta_down *
	      (exp(eta_down*beta1try)/gam_down) *
	      exp(0.5*tau_eps*pow(beta[g][jstar],2)) / sqrt(tau_eps/6.28318);
	  }
	}
	else if(uu<wtc[0]+wtc[2]) {
	  zg_try = 2;
	  if(beta1try<0) acc_prob = 0;
	  else{
	    if(zg[g]==0) acc_prob = (eta_up/eta_down) * 
	      exp(-eta_up*beta1try-eta_down*beta[g][jstar]) * (1./ratiogam) * 
	      pow(beta1try*eta_up,lambda_up-1)/pow(-beta[g][jstar]*eta_down,lambda_down-1);
	    else if(zg[g]==2) acc_prob = pow(beta1try/beta[g][jstar],lambda_up-1) * 
	      exp(-eta_up*(beta1try-beta[g][jstar]));
	    else acc_prob = pow(beta1try*eta_up,lambda_up-1) * eta_up *
	      (exp(-eta_up*beta1try)/gam_up) *
	      exp(0.5*tau_eps*pow(beta[g][jstar],2)) / sqrt(tau_eps/6.28318);
	  }
	}
	else {
	  zg_try = 1;
	  if(zg[g]==0) acc_prob = exp(-0.5*tau_eps*pow(beta1try,2)) *
	    sqrt(tau_eps/6.28318) *
	    exp(-eta_down*beta[g][jstar])*gam_down / 
	    ( eta_down * pow(-beta[g][jstar]*eta_down,lambda_down-1) );
	  else if(zg[g]==2) acc_prob = exp(-0.5*tau_eps*pow(beta1try,2)) *
	    sqrt(tau_eps/6.28318) *
	    exp(eta_up*beta[g][jstar])*gam_up / 
	    ( eta_up * pow(beta[g][jstar]*eta_up,lambda_up-1) );
	  else acc_prob = exp(-0.5*tau_eps*( pow(beta1try,2)-pow(beta[g][jstar],2) ));
	}

//std::cerr << zg[g] << " " << zg_try << " " << acc_prob << "    ";

	uu = rand.Uniform();
	n_try++;
	if(uu<acc_prob){
		nalloc[zg_try]++;
		nalloc[zg[g]]--;
		zg[g] = zg_try;
		beta[g][jstar] = beta1try;
		n_acc++; 
	}

}

//std::cerr << "sum alloc" << nalloc[0]+nalloc[1]+nalloc[2] << " ";

}


/////////////////////////////////////////////////////////////////////////
//// Gibbs for eta_down, eta_up (hyperparam of Uniform in mixture prior)

void update_eta_unif(double &eta_up, double &eta_down, double &aa_eta, double &bb_eta, 
int* nalloc, Random &rand){

  double eta_up_try,eta_down_try,lacc_up,lacc_down,uu;

  eta_up_try = rand.Uniform(aa_eta,bb_eta);
  eta_down_try = rand.Uniform(aa_eta,bb_eta);

  lacc_up = nalloc[2]*(log(eta_up)-log(eta_up_try));
  lacc_down = nalloc[0]*(log(eta_down)-log(eta_down_try));

  uu = rand.Uniform();
  if(uu<exp(lacc_up)) eta_up = eta_up_try;

  uu = rand.Uniform();
  if(uu<exp(lacc_down)) eta_down = eta_down_try;

  //  std::cerr << "  eta " << eta_up << " " << eta_down << std::endl;

}

/////////////////////////////////////////////////////////////////////////
//// Gibbs for eta_down, eta_up (hyperparam of Gamma in mixture prior)

void update_eta(double &eta_up, double &eta_down, double &lambda_up, double &lambda_down, double &aa_eta, double &bb_eta, 
int* zg, double** beta, int* nalloc, int &ngenes, int &nconds, int &jstar, Random &rand){

double dum0,dum2;

//std::cerr << "eta" << " ";

dum0 = 0;
dum2 = 0;
for(int g=0; g<ngenes; ++g) {
	if(zg[g]==0) dum0 += beta[g][jstar];
	if(zg[g]==2) dum2 += beta[g][jstar];
}


eta_down = rand.Gamma( aa_eta + nalloc[0]*lambda_down , bb_eta - dum0 );
eta_up = rand.Gamma( aa_eta + nalloc[2]*lambda_up , bb_eta + dum2 );

// std::cerr << "   dum0 " << dum0 << " nalloc[0] " << nalloc[0];
// std::cerr << "   dum2 " << dum2 << " nalloc[2] " << nalloc[2];
// std::cerr << "    eta's " << eta_down << " " << eta_up << std::endl;

}

/////////////////////////////////////////////////////////////////////////
//// Gibbs for lambda_down, lambda_up (hyperparam of Gamma in mixture prior)

void update_lambda(double &lambda_up, double &lambda_down, double &eta_up, double &eta_down, double &lam1, double &lam2, int &nlam, int* zg, double** beta, int* nalloc, int &ngenes, int &jstar, Random &rand){

  double * lplam0 = new double[nlam];
  double * lplam2 = new double[nlam];
  double * lam_new = new double[nlam];
  double psum0, psum2, uu, pdum;
  int updated;
  
  psum0 = 0;
  psum2 = 0;
  for(int j=0; j<nlam; ++j){
    lam_new[j] = lam1 + j*(lam2-lam1)/float(nlam-1);
    lplam0[j] = 0;
    lplam2[j] = 0;
    for(int g=0; g<ngenes; ++g) {
      if(zg[g]==0) lplam0[j] += (lam_new[j]-1)*log(-beta[g][jstar]);
      if(zg[g]==2) lplam2[j] += (lam_new[j]-1)*log(beta[g][jstar]);
    }
  /////////// replaced gsl with R fns
      // lplam0[j] += nalloc[0]*( lam_new[j]*log(eta_down) - log(gsl_sf_gamma(lam_new[j])) );
      // lplam2[j] += nalloc[2]*( lam_new[j]*log(eta_up) - log(gsl_sf_gamma(lam_new[j])) );
    lplam0[j] += nalloc[0]*( lam_new[j]*log(eta_down) - log(gammafn(lam_new[j])) );
    lplam2[j] += nalloc[2]*( lam_new[j]*log(eta_up) - log(gammafn(lam_new[j])) );
    psum0 += exp(lplam0[j]);
    psum2 += exp(lplam2[j]);
  }

  uu = rand.Uniform()*psum0;
  updated = 0;
  pdum = 0;
  for(int j=0; j<nlam; ++j){
    if(updated==0){
      pdum +=  exp(lplam0[j]);
      if(uu <= pdum){
	lambda_down = lam_new[j];
	updated = 1;
      }
    }
  }
  
  uu = rand.Uniform()*psum2;
  updated = 0;
  pdum = 0;
  for(int j=0; j<nlam; ++j){
    if(updated==0){
      pdum +=  exp(lplam2[j]);
      if(uu <= pdum){
	lambda_up = lam_new[j];
	updated = 1;
      }
    }
  }
  
}

/////////////////////////////////////////////////////////////////////////
//// Gibbs for weights 

void update_wtc(double* wtc, int* nalloc, double* nu, int &ncomps,
Random &rand){

//std::cerr << "wtc" << " ";

// Dirichlet fn overwrites 1st arg
// 1st arg must be type valarray, not usual pointer

std::valarray<double> dumvector(0.0,ncomps);
for(int j=0; j<ncomps; ++j) dumvector[j] = nalloc[j]+nu[j];

rand.Dirichlet(dumvector,ncomps);

for(int j=0; j<ncomps; ++j) wtc[j] = dumvector[j];
///std::cerr << wtc[0] << " " << wtc[1] << " " << wtc[2] << " ";

  // std::cerr << "  wtc " << wtc[1] << std::endl;

}

/////////////////////////////////////////////////////////////////////////
//// Gibbs for tau (Gamma prior)

void update_tau(double** beta, double** tau, double** gamma,
		double** xx, int* indtau, double** ybar, double** ss, double** ydata, 
		double* aa, double* bb, int &like_choice, int &ngenes, int &nconds, 
		int &ntau, int* nreps, int &neffects, Random &rand){

double dum2,dum3;
double sumnreps;
int ind;

//std::cerr << "tau" << " ";

for(int g=0; g<ngenes; ++g) {
for(int t=0; t<ntau; ++t) {
  dum3 = 0;
  sumnreps = 0;
  for(int c=0; c<nconds; ++c) {
    if(indtau[c]==t){			
      dum2 = 0;
      for(int j=0; j<neffects; ++j) {
	dum2 += beta[g][j]*xx[j][c];
      }
      if(like_choice==1){
	dum3 += (nreps[c]-1)*ss[g][c] + nreps[c]*pow(ybar[g][c] - dum2, 2);
      }
      else if(like_choice==2){
	for(int r=0; r<nreps[c]; ++r) {
	  if(c==0) ind = r;
	  else ind = c*nreps[c-1] + r;
	  dum3 += gamma[g][ind]*pow(ydata[g][ind] - dum2, 2);
	}
      }
      sumnreps += nreps[c];
    }
  }
  tau[g][t] = rand.Gamma( aa[t] + sumnreps/2. , bb[t] + dum3/2. );
}}

// std::cerr << "  tau " << tau[0][0] << std::endl;

}

/////////////////////////////////////////////////////////////////////////
//// MH for tau (logNormal prior)

void update_tau_logNorm(double** beta, double** tau, double** gamma,
			double** xx, int* indtau, double** ybar, double** ss, double** ydata, 
			double* aa, double* bb, 
			int &n_tau_lognorm_acc, int &n_tau_lognorm_try,
			int &like_choice, int &ngenes, int &nconds, int &ntau, int* nreps, 
			int &neffects, Random &rand){
  
double dum2,dum3,tau_try,log_acc_prob,uu;
double sumnreps;
int ind;

//std::cerr << "tau" << " ";

for(int g=0; g<ngenes; ++g) {
for(int t=0; t<ntau; ++t) {
  dum3 = 0;
  sumnreps = 0;
  for(int c=0; c<nconds; ++c) {
    if(indtau[c]==t){			
      dum2 = 0;
      for(int j=0; j<neffects; ++j) {
	dum2 += beta[g][j]*xx[j][c];
      }
      if(like_choice==1){
	dum3 += (nreps[c]-1)*ss[g][c] + nreps[c]*pow(ybar[g][c] - dum2, 2);
      }
      else if(like_choice==2){
	for(int r=0; r<nreps[c]; ++r) {
	  if(c==0) ind = r;
	  else ind = c*nreps[c-1] + r;
	  dum3 += gamma[g][ind]*pow(ydata[g][ind] - dum2, 2);
	}
      }
      sumnreps += nreps[c];
    }
  }

  tau_try = rand.Gamma( sumnreps/2. , dum3/2. );
   
  log_acc_prob = ( pow(log(tau[g][t])-aa[t],2) - pow(log(tau_try)-aa[t],2) )*bb[t]/2.;

  uu = rand.Uniform();
  n_tau_lognorm_try++;
  if(uu<exp(log_acc_prob)){
    tau[g][t] = tau_try;
    n_tau_lognorm_acc++; 
  }   
  
}}




// std::cerr << "  tau " << tau[0][0] << std::endl;

}

/////////////////////////////////////////////////////////////////////////
//// Gibbs for tau (with cut between tau and ybar) only with Gamma prior and only Normal likelihood

void update_tau_cut(double** tau, double** ss, int* indtau, double* aa, double* bb, 
int &ngenes, int &nconds, int &ntau, int* nreps, Random &rand){

double dum1,dum3;
double sumnreps;

//std::cerr << "tau" << " ";

for(int g=0; g<ngenes; ++g) {
for(int t=0; t<ntau; ++t) {
   dum3 = 0;
   sumnreps = 0;
   for(int c=0; c<nconds; ++c) {
	if(indtau[c]==t){			
	   dum1 = (nreps[c]-1)*ss[g][c];
	   dum3 = dum3 + dum1;
	   sumnreps += nreps[c];
	}
   }
   tau[g][t] = rand.Gamma( aa[t] + sumnreps/2. , bb[t] + dum3/2. );
}}

}

/////////////////////////////////////////////////////////////////////////
//// Gibbs for b (hyperparam of Gamma) - when gene variances are inverse gamma

void update_bb(double** tau, double* aa, double* bb, double &gg, double &hh, int &ngenes, int &nconds, int &ntau, Random &rand){

double dum1;

//std::cerr << "bb" << " ";

for(int t=0; t<ntau; ++t) {

	dum1 = 0;
	for(int g=0; g<ngenes; ++g) {
		dum1 += tau[g][t];
	}

	bb[t] = rand.Gamma( gg + ngenes*aa[t] , hh + dum1 );

}

// std::cerr << "  bb " << bb[0] << std::endl;

}

/////////////////////////////////////////////////////////////////////////
//// Gibbs for b (hyperparam of LogNorm) - when gene variances are log Normal

void update_bb_logNorm(double** tau, double* aa, double* bb, double &gg, double &hh, 
		       int &ngenes, int &nconds, int &ntau, Random &rand){
  
double dum1;

//std::cerr << "bb" << " ";

for(int t=0; t<ntau; ++t) {

	dum1 = 0;
	for(int g=0; g<ngenes; ++g) {
	  dum1 += pow(log(tau[g][t])-aa[t],2);
	}

	bb[t] = rand.Gamma( gg + ngenes/2. , hh + dum1/2. );

}

// std::cerr << "  bb " << bb[0] << std::endl;

}


/////////////////////////////////////////////////////////////////////////
//// RW MH for a (hyperparam of Gamma) - when gene variances are Inverse Gamma

void update_aa(double &sig_aa, double** tau, double* aa, double* bb, double &gg, double &hh, int &n_aa_acc, int &n_aa_try, int &ngenes, int &nconds, int &ntau, Random &rand){

double logbtau,aa_try,log_acc_prob,uu;

for(int t=0; t<ntau; ++t) {

	/// geom mean of b*tau
	logbtau = 0;
	for(int g=0; g<ngenes; ++g) {
		logbtau += log(tau[g][t]);
	}
	logbtau = logbtau/float(ngenes);
	logbtau += log(bb[t]);

	aa_try = rand.Normal( aa[t],sig_aa );

  /////////// replaced gsl with R fns
	  // log_acc_prob = (gg-1)*(log(aa_try)-log(aa[t])) + 
	  // ngenes*( log(gsl_sf_gamma(aa[t])) - log(gsl_sf_gamma(aa_try)) )
          //     + (aa_try-aa[t])*(ngenes*logbtau-hh);
	log_acc_prob = (gg-1)*(log(aa_try)-log(aa[t])) + 
	  ngenes*( log(gammafn(aa[t])) - log(gammafn(aa_try)) )
               + (aa_try-aa[t])*(ngenes*logbtau-hh);
	
	uu = rand.Uniform();
	n_aa_try++;
	if(uu<exp(log_acc_prob)){
		aa[t] = aa_try;
		n_aa_acc++; 
	}

}

// std::cerr << "  aa " << aa[0] << std::endl;

}

/////////////////////////////////////////////////////////////////////////
//// Gibbs for a (hyperparam of logNormal) - when gene variances are logNormal

void update_aa_logNorm(double &sig_aa, double** tau, double* aa, double* bb, double &tau_var, 
		      int &ngenes, int &nconds, int &ntau, Random &rand){

  double logbtau;

 for(int t=0; t<ntau; ++t) {

   /// geom mean of b*tau
   logbtau = 0;
   for(int g=0; g<ngenes; ++g) {
     logbtau += log(tau[g][t]);
   }
   
   aa[t] = rand.Normal( bb[t]*logbtau/( ngenes*bb[t] + tau_var ) , 1.0/sqrt(ngenes*bb[t] + tau_var) );
 }

// std::cerr << "  aa " << aa[0] << std::endl;

}


/////////////////////////////////////////////////////////////////////////
//// Gibbs for gamma (scale param for t as mix of Normals)

void update_gamma(double** beta, double** tau, double** gamma,
		  double** xx, int* indtau, double** ydata, double* df, 
		  int &ngenes, int &nconds, int &ntau, int* nreps, 
		  int &neffects, Random &rand){

double dum2,dum3;
int ind;

for(int g=0; g<ngenes; ++g) {
for(int c=0; c<nconds; ++c) {
  dum2 = 0;
  for(int j=0; j<neffects; ++j) {
    dum2 += beta[g][j]*xx[j][c];
  }
  for(int r=0; r<nreps[c]; ++r) {
    if(c==0) ind = r;
    else ind = c*nreps[c-1] + r;
    dum3 = pow(ydata[g][ind] - dum2, 2);
    gamma[g][ind] = rand.Gamma( (df[c]+1)/2. , df[c]/2. + tau[g][indtau[c]]*dum3/2. );
  }
}}

}

/////////////////////////////////////////////////////////////////////////
//// mixed/posterior predictive pvals for ss_pred and ybar_pred (NORMAL LIKELIHOOD ONLY)

void predict(double** ybar_pred1, double** ybar_pred2, double** ybar_pred3, double** ybar_pred4, 
double** ss_pred1, double** ss_pred2, double** pval_post_ss, double** pval_mix_ss, 
double** pval_post_ybar, double** pval_mix1_ybar, double** pval_mix2_ybar, double** pval_mix3_ybar, 
double** pval_partial_ss, double* pval_partial_ybar, double** norm_ss, double* norm_ybar,
double** ybar, double** ss, 
double** tau, double** gamma, double* aa, double* bb, 
int* zg, double** beta, double** xx, int* indtau, double* wtc, double &tau_eps,
double &eta_up, double &eta_down, double &lambda_up, double &lambda_down, int &jstar, 
int &move_choice_bz, int &move_choice_tau, int &ngenes, 
int &nconds, int* nreps, int &neffects, Random &rand){

double tau_pred,dum2,dum2_pred,log_tau_pred;
double beta_jstar_pred=0;
double taunew,den_ss,den_ybar;

for(int g=0; g<ngenes; ++g) {
if(jstar!=-1){
	/// predict new beta from mixture
	//uu = rand.Uniform();
	if(move_choice_bz==1 | move_choice_bz==2){
		//if(uu<wtc[0]) {
		if(zg[g]==0){
			beta_jstar_pred = rand.Uniform(-eta_down,0);
		}
		//else if(uu<wtc[0]+wtc[2]) {
		else if(zg[g]==2){
			beta_jstar_pred = rand.Uniform(0,eta_up);
		}
		else {
			beta_jstar_pred = 0;
		}
	}
	else if(move_choice_bz==3 | move_choice_bz==4 |  move_choice_bz==5){
		//if(uu<wtc[0]) {
		if(zg[g]==0){
			beta_jstar_pred = -rand.Gamma( lambda_down, eta_down );
		}
		//else if(uu<wtc[0]+wtc[2]) {
		else if(zg[g]==2){
			beta_jstar_pred = rand.Gamma( lambda_up, eta_up );
		}
		else {
			if(move_choice_bz==3 | move_choice_bz==4) beta_jstar_pred = 0;
                        else if(move_choice_bz==5) beta_jstar_pred = rand.Normal(0,1./sqrt(tau_eps));
		}
	}
}
for(int c=0; c<nconds; ++c) {

	/// get dum2
	dum2 = 0;
	for(int jp=0; jp<neffects; ++jp) {
		if(jp!=jstar) dum2 += beta[g][jp]*xx[jp][c];
	}
	if(jstar!=-1){
		dum2_pred = dum2 + beta_jstar_pred*xx[jstar][c];
		dum2 += beta[g][jstar]*xx[jstar][c];
	}
	else dum2_pred = dum2;

	/// predict new tau
	if(move_choice_tau==1) tau_pred = rand.Gamma( aa[indtau[c]] ,  bb[indtau[c]] );
	else if(move_choice_tau==2){
	  log_tau_pred = rand.Normal( aa[indtau[c]] ,  1.0/sqrt(bb[indtau[c]]) );
	  tau_pred = exp(log_tau_pred);
	}
	else{
	  tau_pred = 0;
	  std::cerr << "  move choice invalid for tau " << std::endl;
	}

	/// posterior predictive
	ybar_pred1[g][c] = rand.Normal( dum2, 
             sqrt(1.0/(nreps[c]*tau[g][indtau[c]])) );
	ss_pred1[g][c] = rand.Gamma( (nreps[c]-1)/2.0 , 
             (nreps[c]-1)*tau[g][indtau[c]]/2.0 ); 

	/// mixed predictive (predicting tau only)
	ybar_pred2[g][c] = rand.Normal( dum2, 
                sqrt(1.0/(nreps[c]*tau_pred)) );
	ss_pred2[g][c] = rand.Gamma( (nreps[c]-1)/2.0 , 
               (nreps[c]-1)*tau_pred/2.0 ); 

	/// mixed predictive (predicting beta only)
	ybar_pred3[g][c] = rand.Normal( dum2_pred, 
                sqrt(1.0/(nreps[c]*tau[g][indtau[c]])) );

	/// mixed predictive (predicting both tau and beta)
	ybar_pred4[g][c] = rand.Normal( dum2_pred, 
                sqrt(1.0/(nreps[c]*tau_pred)) );

	if(ss_pred1[g][c] >= ss[g][c]) pval_post_ss[g][c] += 1;
	if(ss_pred2[g][c] >= ss[g][c]) pval_mix_ss[g][c] += 1;

	/// partial pvals for S2
	  den_ss = pow(tau[g][indtau[c]],0.5*(nreps[c]-1)) * pow(ss[g][c],0.5*(nreps[c]-3)) * 
	    exp(-0.5*(nreps[c]-1)*ss[g][c]*tau[g][indtau[c]]);
	  norm_ss[g][c] += 1.0/den_ss;
	  if(ss_pred1[g][c] > ss[g][c]) pval_partial_ss[g][c] += 1.0/den_ss;

}
/// !!! only valid for alpha,delta parametrization (not for paired data)

/// delta conditioned on zg
if(ybar_pred1[g][1]-ybar_pred1[g][0] >= ybar[g][1]-ybar[g][0]) pval_post_ybar[g][zg[g]] += 1;
if(ybar_pred2[g][1]-ybar_pred2[g][0] >= ybar[g][1]-ybar[g][0]) pval_mix1_ybar[g][zg[g]] += 1;
if(ybar_pred3[g][1]-ybar_pred3[g][0] >= ybar[g][1]-ybar[g][0]) pval_mix2_ybar[g][zg[g]] += 1;
if(ybar_pred4[g][1]-ybar_pred4[g][0] >= ybar[g][1]-ybar[g][0]) pval_mix3_ybar[g][zg[g]] += 1;

/// partial p-vals for diff
  taunew = nreps[0]*nreps[1]*tau[g][indtau[0]]*tau[g][indtau[1]] / 
    (nreps[0]*tau[g][indtau[0]] + nreps[1]*tau[g][indtau[1]]);
  den_ybar = sqrt(taunew)*exp(-0.5*taunew*pow(ybar[g][1]-ybar[g][0]-beta[g][jstar],2));
  norm_ybar[g] += 1.0/den_ybar;
  if(ybar_pred1[g][1]-ybar_pred1[g][0] > ybar[g][1]-ybar[g][0]) pval_partial_ybar[g] += 1.0/den_ybar;

} // end of g-loop

}

/////////////////////////////////////////////////////////////////////////
//// mixed/posterior predictive pvals for ss_pred and ybar_pred (FOR T LIKELIHOOD)

void predict_t(double** ybar_pred1, double** ybar_pred2, double** ybar_pred3, double** ybar_pred4, 
double** ss_pred1, double** ss_pred2, double** pval_post_ss, double** pval_mix_ss, 
double** pval_post_ybar, double** pval_mix1_ybar, double** pval_mix2_ybar, double** pval_mix3_ybar, 
double** ybar, double** ss, double** tau, double** gamma, double* df, double* aa, double* bb, 
int* zg, double** beta, double** xx, int* indtau, double* wtc, double &tau_eps,
double &eta_up, double &eta_down, double &lambda_up, double &lambda_down, int &jstar, 
int &move_choice_bz, int &move_choice_tau, int &ngenes, 
int &nconds, int* nreps, int &neffects, Random &rand){

double tau_pred,dum2,dum2_pred,log_tau_pred,gamma_pred;
double beta_jstar_pred=0;
double ydata_pred1,ydata_pred2,ydata_pred3,ydata_pred4;
int ind;

for(int g=0; g<ngenes; ++g) {
if(jstar!=-1){
	/// predict new beta from mixture
	//uu = rand.Uniform();
	if(move_choice_bz==1 | move_choice_bz==2){
		//if(uu<wtc[0]) {
		if(zg[g]==0){
			beta_jstar_pred = rand.Uniform(-eta_down,0);
		}
		//else if(uu<wtc[0]+wtc[2]) {
		else if(zg[g]==2){
			beta_jstar_pred = rand.Uniform(0,eta_up);
		}
		else {
			beta_jstar_pred = 0;
		}
	}
	else if(move_choice_bz==3 | move_choice_bz==4 |  move_choice_bz==5){
		//if(uu<wtc[0]) {
		if(zg[g]==0){
			beta_jstar_pred = -rand.Gamma( lambda_down, eta_down );
		}
		//else if(uu<wtc[0]+wtc[2]) {
		else if(zg[g]==2){
			beta_jstar_pred = rand.Gamma( lambda_up, eta_up );
		}
		else {
			if(move_choice_bz==3 | move_choice_bz==4) beta_jstar_pred = 0;
                        else if(move_choice_bz==5) beta_jstar_pred = rand.Normal(0,1./sqrt(tau_eps));
		}
	}
}
for(int c=0; c<nconds; ++c) {

	/// get dum2
	dum2 = 0;
	for(int jp=0; jp<neffects; ++jp) {
		if(jp!=jstar) dum2 += beta[g][jp]*xx[jp][c];
	}
	if(jstar!=-1){
		dum2_pred = dum2 + beta_jstar_pred*xx[jstar][c];
		dum2 += beta[g][jstar]*xx[jstar][c];
	}
	else dum2_pred = dum2;

	/// predict new tau
	if(move_choice_tau==1) tau_pred = rand.Gamma( aa[indtau[c]] ,  bb[indtau[c]] );
	else if(move_choice_tau==2){
	  log_tau_pred = rand.Normal( aa[indtau[c]] ,  1.0/sqrt(bb[indtau[c]]) );
	  tau_pred = exp(log_tau_pred);
	}
	else{
	  tau_pred = 0;
	  std::cerr << "  move choice invalid for tau " << std::endl;
	}

	/// predict new gamma and ydata
	ybar_pred1[g][c] = 0;
	ybar_pred2[g][c] = 0;
	ybar_pred3[g][c] = 0;
	ybar_pred4[g][c] = 0;
	ss_pred1[g][c] = 0;
	ss_pred2[g][c] = 0;
	for(int r=0; r<nreps[c]; ++r) {
	  if(c==0) ind = r;
	  else ind = c*nreps[c-1] + r;
	  gamma_pred = rand.Gamma( df[c]/2. , df[c]/2. );
	  ydata_pred1 = rand.Normal( dum2, sqrt(1.0/(tau[g][indtau[c]]*gamma_pred)) );
	  ydata_pred2 = rand.Normal( dum2, sqrt(1.0/(tau_pred*gamma_pred)) );
	  ydata_pred3 = rand.Normal( dum2_pred, sqrt(1.0/(tau[g][indtau[c]]*gamma_pred)) );
	  ydata_pred4 = rand.Normal( dum2_pred, sqrt(1.0/(tau_pred*gamma_pred)) );
	  ybar_pred1[g][c] += ydata_pred1;
	  ybar_pred2[g][c] += ydata_pred2;
	  ybar_pred3[g][c] += ydata_pred3;
	  ybar_pred4[g][c] += ydata_pred4;
	  if(r>0){
	    ss_pred1[g][c] = ss_pred1[g][c]*(r-1)/float(r) + 
	      pow(( ydata_pred1 - ybar_pred1[g][c]/float(r+1) )/float(r) ,2)*(r+1);
	    ss_pred2[g][c] = ss_pred2[g][c]*(r-1)/float(r) + 
	      pow(( ydata_pred2 - ybar_pred2[g][c]/float(r+1) )/float(r) ,2)*(r+1);
	  }
	}
	ybar_pred1[g][c] = ybar_pred1[g][c]/float(nreps[c]);
	ybar_pred2[g][c] = ybar_pred2[g][c]/float(nreps[c]);
	ybar_pred3[g][c] = ybar_pred3[g][c]/float(nreps[c]);
	ybar_pred4[g][c] = ybar_pred4[g][c]/float(nreps[c]);
       
	if(ss_pred1[g][c] >= ss[g][c]) pval_post_ss[g][c] += 1;
	if(ss_pred2[g][c] >= ss[g][c]) pval_mix_ss[g][c] += 1;

} // end of c-loop

/// !!! only valid for alpha,delta parametrization (not for paired data)

/// delta conditioned on zg
if(ybar_pred1[g][1]-ybar_pred1[g][0] >= ybar[g][1]-ybar[g][0]) pval_post_ybar[g][zg[g]] += 1;
if(ybar_pred2[g][1]-ybar_pred2[g][0] >= ybar[g][1]-ybar[g][0]) pval_mix1_ybar[g][zg[g]] += 1;
if(ybar_pred3[g][1]-ybar_pred3[g][0] >= ybar[g][1]-ybar[g][0]) pval_mix2_ybar[g][zg[g]] += 1;
if(ybar_pred4[g][1]-ybar_pred4[g][0] >= ybar[g][1]-ybar[g][0]) pval_mix3_ybar[g][zg[g]] += 1;

} // end of g-loop

}


/////////////////////////////////////////////////////////////////////////
//// deviance

void deviance_calc(double &deviance1, double &deviance2, double** beta, double** tau, 
		   double** gamma, double** xx, int* indtau, double** ybar, 
		   double** ss, double** ydata, double* df, int &like_choice, int &ngenes, 
		   int &nconds, int* nreps, int &neffects){

double dum1, dum2;
int ind;

deviance1=0;
deviance2=0;
for(int g=0; g<ngenes; ++g) {
for(int c=0; c<nconds; ++c) {
  dum2 = 0;
  for(int j=0; j<neffects; ++j) {
    dum2 += beta[g][j]*xx[j][c];
  }
  if(like_choice==1){
    dum1 = (nreps[c]-1)*ss[g][c] + nreps[c]*pow(ybar[g][c] - dum2, 2);
    deviance1 += tau[g][indtau[c]]*dum1 - nreps[c]*log(tau[g][indtau[c]]);
  }
  else if(like_choice==2){
    for(int r=0; r<nreps[c]; ++r) {
      if(c==0) ind = r;
      else ind = c*nreps[c-1] + r;
      dum1 = pow(ydata[g][ind] - dum2, 2);
      deviance1 += tau[g][indtau[c]]*gamma[g][ind]*dum1 - log(gamma[g][ind]*tau[g][indtau[c]]);
      deviance2 += (df[c]+1)*log(1 + dum1*tau[g][indtau[c]]/df[c]) - log(df[c]*tau[g][indtau[c]]);
    }
  }

}}


}




