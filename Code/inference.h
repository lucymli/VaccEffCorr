#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <string>

class Data {
  Rcpp::NumericMatrix swab_data_v; // nrow=total number of vaccinated, ncol=swabs
  Rcpp::NumericMatrix swab_data_nv; // nrow=total number of vaccinated, ncol=swabs
  Rcpp::NumericVector swab_times; // timing of swabs
  Rcpp::NumericMatrix ab_data; // nrow=total number of vaccinated, ncol=serotypes in a vaccine
  int n_vacc;
  int n_nvacc;
  int n_swabs;
  int n_vtypes;
  int n_nvtypes;
  int n_tot;
public:
  Data (Rcpp::List, Rcpp::NumericVector, int, int);
};

Data::Data(Rcpp::List input_data, Rcpp::NumericVector timings, int vtypes, int nvtypes) {
  swab_times = timings;
  n_vtypes = vtypes;
  n_nvtypes = nvtypes;
  n_tot = n_vtypes + n_nvtypes;
  swab_data_v = input_data["vdata"];
  swab_data_nv = input_data["nvdata"];
  ab_data = input_data["abdata"];
}

class Param {
  int n_vtypes;
  int n_nvtypes;
  int n_tot;
  Rcpp::NumericVector tempparam;
  Rcpp::NumericVector lambda;
  Rcpp::NumericVector mu;
  Rcpp::NumericVector thetaSI;
  Rcpp::NumericVector thetaIS;
  Rcpp::NumericVector p0;
  double frailtySI;
  double frailtyIS;
  double interaction;
  int nblocks;
  int block_ptr;
  Rcpp::NumericVector acceptance;
  Rcpp::NumericVector variance;
  double llik;
  double lprior;
public:
  Param (Rcpp::List);
  void next_block ();
  double calc_llik (Data);
  double calc_lprior ();
  double uni_propose (double, double, int) const;
  void alter_param (bool);
  void mcmc_move (Data);
  void print_params (double *) const;
};

Param::Param (Rcpp::List input) {
  n_vtypes = input["ntypes"];
  n_tot = input["ntot"];
  n_nvtypes = n_tot-n_vtypes;
  tempparam = Rcpp::NumericVector(1);
  lambda = Rcpp::NumericVector(input["lambda"]);
  mu = Rcpp::NumericVector(input["mu"]);
  thetaSI = Rcpp::NumericVector(input["thetaSI"]);
  thetaIS = Rcpp::NumericVector(input["thetaIS"]);
  p0 = Rcpp::NumericVector(input["p0"]);
  frailtySI = input["frailtySI"];
  frailtyIS = input["frailtyIS"];
  interaction = input["interaction"];
  nblocks = 8;
  block_ptr = 0;
  acceptance = Rcpp::NumericVector(nblocks);
  variance = input["variances"];
  llik = 0.0;
  lprior = Rcpp::NumericVector(lprior);
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
        for (int i=0; i<n_tot; i++) {
          mu[i] = uni_propose(tempparam[i], variance[block_ptr], ntries);
        }
      }
    }
    case 3: {
      if (reject) thetaSI = tempparam; return;
      tempparam = thetaSI;
      for (int i=0; i<n_vtypes; i++) {
        for (int i=0; i<n_tot; i++) {
          thetaSI[i] = uni_propose(tempparam[i], variance[block_ptr], ntries);
        }
      }
    }
    case 4: {
      if (reject) thetaIS = tempparam; return;
      tempparam = thetaIS;
      for (int i=0; i<n_vtypes; i++) {
        for (int i=0; i<n_tot; i++) {
          thetaIS[i] = uni_propose(tempparam[i], variance[block_ptr], ntries);
        }
      }
    }
    case 5: {
      if (reject) p0 = tempparam; return;
      tempparam = p0;
      for (int i=0; i<n_vtypes; i++) {
        for (int i=0; i<n_tot; i++) {
          p0[i] = uni_propose(tempparam[i], variance[block_ptr], ntries);
        }
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

void Param::mcmc_move (Data data) {
  // Code to propose
  alter_param (false);
  double newllik = calc_llik(data);
  double newlprior = calc_lprior();
  double z = log(R::unif_rand());
  if (z > (newllik+newlprior-llik-lprior)) alter_param (true);
  next_block();
}


void Param::print_params (double * output) const {
}