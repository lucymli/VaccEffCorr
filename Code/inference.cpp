#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <string>
#include <stdio.h>      /* printf */
#include <math.h>       /* log */
using namespace Rcpp;

class Param {
  int n_vtypes, n_nvtypes, n_tot, n_blocks, block_ptr; // counts
  std::vector <double> tempparam, lambda, mu; // parameter vectors
  std::vector <double> thetaSI, thetaIS, p0;
  double frailtySI, frailtyIS, interaction; // non-serotype-specific parameters
  std::vector <double> accepted, rejected, variance; // mcmc vectors
  double llik; // overall log likelihood
  std::vector <double> llik_vec; // log likelihood for each block
  double lprior; // overall prior
  std::vector <double> lprior_vec; // prior for each block
  arma::mat transitions;
public:
  Param (int, int, std::vector<double>, std::vector<double>);
  void next_block ();
  double calc_llik (Data, int) const;
  double calc_llik (Data) const;
  double calc_lprior (int) const;
  double calc_lprior () const;
  double uni_propose (double, double, int) const;
  void alter_param (bool);
  void mcmc_move (Data, bool, double);
  double operator [] (int);
  void print_params (std::vector <double> &);
};


Param::Param (int ntypes, int ntot, std::vector <double> inputs, std::vector <double> input_var) {
  n_vtypes = ntypes;
  n_tot = ntot;
  n_nvtypes = n_tot-n_vtypes;
  tempparam = std::vector <double> (n_tot);
  frailtySI = inputs[n_tot*2+n_vtypes*3];
  frailtyIS = inputs[n_tot*2+n_vtypes*3+1];
  interaction = inputs[n_tot*2+n_vtypes*3+2];
  for (int i=0; i<n_vtypes; i++) {
    lambda.push_back(inputs[i]);
    mu.push_back(inputs[i+n_tot]);
    thetaSI.push_back(inputs[i+n_tot*2]);
    thetaIS.push_back(inputs[i+n_tot*2+n_vtypes]);
    p0.push_back(inputs[i+n_tot*2+n_vtypes*2]);
  }
  for (int i=n_vtypes; i<n_tot; i++) {
    lambda.push_back(inputs[i]);
    mu.push_back(inputs[i+n_tot]);
  }
  n_blocks = 8;
  block_ptr = 0;
  accepted = std::vector <double>(n_blocks);
  rejected = std::vector <double>(n_blocks);
  variance = input_var;
  llik = 0.0;
  lprior = 0.0;
  arma::mat;
}

void Param::next_block() {
  block_ptr++;
  if (block_ptr >= n_blocks) block_ptr = 0;
}

double Param::calc_llik (Data data, int lik_type) {
  double newllik = 0.0;
  return (newllik);
}

double Param::calc_llik (Data data) const {
  double newllik = 0.0;
  newllik += calc_llik(data, 0);
  return (newllik);
}

double Param::calc_lprior (int block_i) const {
  double newlprior = 0.0;
  switch (block_i) {
  case 0: for (auto i=lambda.begin(); i<lambda.end(); i++) newlprior += R::dlnorm(*i, -5.9, 0.2, 1);
  case 1: for (auto i=mu.begin(); i<mu.end(); i++) newlprior += R::dlnorm(*i, -5.9, 0.2, 1);
  case 2: for (auto i=thetaSI.begin(); i<thetaSI.end(); i++) newlprior += R::dunif(*i, 0, 1, 1);
  case 3: for (auto i=thetaIS.begin(); i<thetaIS.end(); i++) newlprior += R::dunif(*i, 0.5, 1, 1);
  case 4: for (auto i=p0.begin(); i<p0.end(); i++) newlprior += R::dunif(*i, 0, 1, 1);
  case 5: newlprior += R::dexp(frailtySI, 1, 1);
  case 6: newlprior += R::dexp(frailtyIS, 1, 1);
  case 7: newlprior += R::dunif(interaction, 0, 1, 1);
  }
  return (newlprior);
}

double Param::calc_lprior () const {
  double newlprior = 0.0;
  for (int i=0; i<n_blocks; i++) {
    newlprior += calc_lprior(i);
  }
  return (newlprior);
}

double Param::uni_propose (double oldval, double var, int ntries) const {
  double newval = R::rnorm(oldval, var);
  while (newval < 0) {
    newval = R::rnorm(oldval, var);
    ntries--;
    if (ntries < 0) break;
  }
  return (newval);
}

void Param::alter_param (bool reject) {
  int ntries = 10;
  switch (block_ptr) {
  case 1: {
    if (reject) lambda = tempparam; return;
    tempparam = lambda;
    for (int i=0; i<n_tot; i++) {
      lambda[i] = uni_propose(tempparam[i], variance[block_ptr], ntries);
    }
  }
  case 2: {
    if (reject) mu = tempparam; return;
    tempparam = mu;
    for (int i=0; i<n_tot; i++) {
      mu[i] = uni_propose(tempparam[i], variance[block_ptr], ntries);
    }
  }
  case 3: {
    if (reject) thetaSI = tempparam; return;
    tempparam = thetaSI;
    for (int i=0; i<n_vtypes; i++) {
      thetaSI[i] = uni_propose(tempparam[i], variance[block_ptr], ntries);
    }
  }
  case 4: {
    if (reject) thetaIS = tempparam; return;
    tempparam = thetaIS;
    for (int i=0; i<n_vtypes; i++) {
      thetaIS[i] = uni_propose(tempparam[i], variance[block_ptr], ntries);
    }
  }
  case 5: {
    if (reject) p0 = tempparam; return;
    tempparam = p0;
    for (int i=0; i<n_vtypes; i++) {
      p0[i] = uni_propose(tempparam[i], variance[block_ptr], ntries);
    }
  }
  case 6: {
    if (reject) frailtySI = tempparam[0]; return;
    tempparam[0] = frailtySI;
    frailtySI = uni_propose(tempparam[0], variance[block_ptr], ntries);
  }
  case 7: {
    if (reject) frailtyIS = tempparam[0]; return;
    tempparam[0] = frailtyIS;
    frailtyIS = uni_propose(tempparam[0], variance[block_ptr], ntries);
  }
  case 8: {
    if (reject) interaction = tempparam[0]; return;
    tempparam[0] = interaction;
    interaction = uni_propose(tempparam[0], variance[block_ptr], ntries);
  }
  }
}

void Param::mcmc_move (Data data, bool adapt, double optimal_adapt) {
  // Code to propose
  alter_param (false);
  double newllik = calc_llik(data);
  double newlprior = calc_lprior();
  double z = R::unif_rand();
  bool reject = log(z) > (newllik+newlprior-llik-lprior);
  if (reject) {
    alter_param (true);
    rejected[block_ptr]++;
  } else {
    accepted[block_ptr]++;
  }
  if (adapt) {
    double change = exp(0.999/2.0*(pow(variance[block_ptr], 0.5)-optimal_adapt));
    variance[block_ptr] *= change*change;
    std::fill(accepted.begin(), accepted.end(), 0.0);
    std::fill(rejected.begin(), rejected.end(), 0.0);
  }
  next_block();
}


double Param::operator[](int i) {
  double value;
  if (i < n_tot) value = lambda[i];
  else if (i < n_tot*2) value = mu[i-n_tot*2];
  else if (i < n_tot*2+n_vtypes) value = thetaSI[i-n_tot*2-n_vtypes];
  else if (i < n_tot*2+n_vtypes*2) value = thetaIS[i-n_tot*2-n_vtypes*2];
  else if (i < n_tot*2+n_vtypes*3) value = p0[i-n_tot*2-n_vtypes*3];
  else if (i == n_tot*2+n_vtypes*3) value = frailtySI;
  else if (i == n_tot*2+n_vtypes*3+1) value = frailtyIS;
  else value = interaction;
  return (value);
}

void Param::print_params (std::vector <double> & output) {
  for (int i=0; i<(n_tot*2+n_vtypes*3+3); i++) {
    output[i] = (*this)[i];
  }
}


class Data {
  std::vector <double> swab_data_v; // nrow=total number of vaccinated, ncol=swabs
  std::vector <double> swab_data_nv; // nrow=total number of vaccinated, ncol=swabs
  std::vector <double> swab_times; // timing of swabs
  std::vector <double> ab_data; // nrow=total number of vaccinated, ncol=serotypes in a vaccine
  int n_vacc;
  int n_nvacc;
  int n_swabs;
  int n_vtypes;
  int n_nvtypes;
  int n_tot;
public:
  Data (std::vector<double>, std::vector<double>, std::vector<double>,
        std::vector <double>, int, int);
};


Data::Data (std::vector<double> vdata, std::vector<double>nvdata, std::vector<double>ab_data,
            std::vector <double>timings, int vtypes, int nvtypes) {
  swab_times = timings;
  n_vtypes = vtypes;
  n_nvtypes = nvtypes;
  n_tot = n_vtypes + n_nvtypes;
  swab_data_v = vdata;
  swab_data_nv = nvdata;
  ab_data = ab_data;
}

// [[Rcpp::export]]
NumericVector mcmc(std::vector <double> params, // vector of parameters
                   std::vector <double> swab_data_v, // nrow=total number of vaccinated, ncol=swabs
                   std::vector <double> swab_data_nv, // nrow=total number of vaccinated, ncol=swabs
                   std::vector <double> swab_times, // timing of swabs
                   std::vector <double> ab_data, // nrow=total number of vaccinated, ncol=serotypes in a vaccine
                   std::vector <double> mcmc_options // options for MCMC algorithm
                   ) {
  int ntypes = vt_parameters.ncol();
  int ntot = nvt_parameters.ncol() + ntypes;
  int n_vt_params = vt_parameters.ncol()*vt_parameters.nrow();
  int n_nvt_params = nvt_parameters.ncol()*nvt_parameters.nrow();
  int n_other_params = params.size();
  int nparams = n_vt_params + n_nvt_params + n_other_params;
  NumericMatrix rates(ntot+1, ntot+1);
  int niter = mcmc_options["niter"];
  int sample_every = mcmc_options["sample_every"];
  std::string filename = mcmc_options["filename"];
  NumericMatrix old_vt_parameters = vt_parameters;
  NumericMatrix old_nvt_parameters = nvt_parameters;
  NumericVector old_params = params;
  NumericMatrix results(niter, nparams+3);
  // Calculate likelihood and prior of initial parameters
  double oldloglik = get_likelihood (vt_parameters, nvt_parameters, params, 
                                     swab_data_v, swab_data_nv, swab_times);
  double oldprior = get_prior(vt_parameters, nvt_parameters, params);
  double newloglik;
  double newprior;
  for (int iter=1; iter<niter; iter++) {
    // 
    // Save parameter values
    if ((iter%sample_every)==0) {
      for (int i=0; i<n_vt_params; i++) {
        results[iter, i] = vt_parameters[i];
      }
      for (int i=0; i<n_nvt_params; i++) {
        results[iter, i+n_vt_params] = nvt_parameters[i];
      }
      for (int i=0; i<n_other_params; ++i) {
        results[iter, i+n_vt_params+n_nvt_params] = params[i];
      }
    }
  }
  // arma::mat rates = 
  return results;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
load(sim.data.RData)
mcmc.options <- list(niter=1000, sample_every=10, filename="testmcmc.txt")
mcmc(42)
*/
