#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <string>
#include <stdio.h>      /* printf */
#include <math.h>       /* log */
#include <vector>
#include "Data.h"

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