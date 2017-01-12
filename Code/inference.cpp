#include "Data.h"
#include "Param.h"
using namespace Rcpp;



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
