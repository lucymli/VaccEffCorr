#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <string>

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
  NumericVector tempparam;
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
  void mcmc_move ();
};

Param::Param (List input) {
  tempparam = NumericVector(1);
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
  
}

void Param::mcmc_move (Data data) {
  // Code to propose
  double newllik = calc_llik(data);
  double newlprior = calc_lprior();
}
