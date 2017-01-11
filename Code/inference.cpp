#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <string>
using namespace Rcpp;

class Data {
  NumericMatrix swab_data_v; // nrow=total number of vaccinated, ncol=swabs
  NumericMatrix swab_data_nv; // nrow=total number of vaccinated, ncol=swabs
  NumericVector swab_times; // timing of swabs
  NumericMatrix ab_data; // nrow=total number of vaccinated, ncol=serotypes in a vaccine
  int n_vacc;
  int n_nvacc;
  int n_swabs;
  int n_vtypes;
  int n_nvtypes;
  int n_tot;
public:
  Data (List, NumericVector, int, int);
};

Data::Data(List input_data, NumericVector timings, int vtypes, int nvtypes) {
  swab_times = timings;
  n_vtypes = vtypes;
  n_nvtypes = nvtypes;
  n_tot = n_vtypes + n_nvtypes;
  swab_data_v = input_data["vdata"];
  swab_data_nv = input_data["nvdata"];
  ab_data = input_data["abdata"];
}

class Param {
  NumericVector lambda;
  NumericVector mu;
  NumericVector thetaSI;
  NumericVector thetaIS;
  NumericVector p0;
  double frailtySI;
  double frailtyIS;
  double interaction;
  int nblocks;
  int block_ptr;
  NumericVector acceptance;
  double llik;
  double lprior;
public:
  Param (List);
  void next_block ();
  double calc_llik (Data);
  double calc_lprior ();
  void propose ();
};

Param::Param (List input) {
  lambda = NumericVector(input["lambda"]);
  mu = NumericVector(input["mu"]);
  thetaSI = NumericVector(input["thetaSI"]);
  thetaIS = NumericVector(input["thetaIS"]);
  p0 = NumericVector(input["p0"]);
  frailtySI = input["frailtySI"];
  frailtyIS = input["frailtyIS"];
  interaction = input["interaction"];
  nblocks = 8;
  block_ptr = 0;
  acceptance = NumericVector(nblocks);
  llik = 0.0;
  lprior = NumericVector(lprior);
}

void Param::next_block() {
  block_ptr++;
  if (block_ptr >= nblocks) block_ptr = 0;
}

double Param::calc_llik (Data data) {
  double newllik;
  return (newllik);
}

double Param::calc_lprior () {
  double newlprior;
  return (newlprior);
}

void Param::propose () {
  // Code to propose
}


// [[Rcpp::export]]
NumericVector mcmc(NumericMatrix vt_parameters, // matrix of parameters specific to each serotype in the vaccine
                   NumericMatrix nvt_parameters, // matrix of parameters specific to each serotype not in the vaccine
                   NumericVector params, // other parameters not specific to serotypes
                   NumericMatrix swab_data_v, // nrow=total number of vaccinated, ncol=swabs
                   NumericMatrix swab_data_nv, // nrow=total number of vaccinated, ncol=swabs
                   NumericVector swab_times, // timing of swabs
                   NumericMatrix ab_data, // nrow=total number of vaccinated, ncol=serotypes in a vaccine
                   List mcmc_options // options for MCMC algorithm
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
