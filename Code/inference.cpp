#include <RcppArmadillo.h>
#include <stdio.h>      /* printf */
#include <math.h>       /* log */
#include <iostream>
#include <fstream>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

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
  int operator[] (std::string) const;
  double get_swab_time (int);
  double get_swab_v (int, int);
  double get_swab_nv (int, int);
};


Data::Data (std::vector<double> vdata, std::vector<double>nvdata, std::vector<double>ab_data,
            std::vector <double>timings, int vtypes, int nvtypes) {
  for (int i=1; i<timings.end()-timings.begin(); i++) swab_times[i] = timings[i]-timings[i-1];
  swab_times[0] = timings[0];
  n_vtypes = vtypes;
  n_nvtypes = nvtypes;
  n_tot = n_vtypes + n_nvtypes;
  swab_data_v = vdata;
  swab_data_nv = nvdata;
  ab_data = ab_data;
}

int Data::operator[] (std::string name) const {
  if (name=="n_vacc") return (n_vacc);
  else if (name=="n_nvacc") return (n_nvacc);
  else if (name=="n_swabs") return (n_swabs);
  else if (name=="n_vtypes") return (n_vtypes);
  else if (name=="n_nvtypes") return (n_nvtypes);
  else return(n_tot);
}

double Data::get_swab_time (int time) {
  return (swab_times[time]);
}

double Data::get_swab_v (int ind_i, int swab_i) {
  return (swab_data_v[swab_i*n_vacc+ind_i]);
}
double Data::get_swab_nv (int ind_i, int swab_i) {
  return (swab_data_nv[swab_i*n_nvacc+ind_i]);
}

class Param {
  int n_vtypes, n_nvtypes, n_tot, n_blocks, block_ptr; // counts
  std::vector <double> tempparam, lambda, mu; // parameter vectors
  std::vector <double> thetaSI, thetaIS, p0;
  double frailtySI, frailtyIS, interaction; // non-serotype-specific parameters
  std::vector <double> accepted, rejected, proposal_sd; // mcmc vectors
  double llik; // overall log likelihood
  std::vector <double> llik_vec; // log likelihood for each block
  double lprior; // overall prior
  std::vector <double> lprior_vec; // prior for each block
  arma::mat transitions;
  arma::mat stationary_prev;
  arma::mat transitions_t;
  std::vector <double> ind_frailty_SI;
  std::vector <double> ind_frailty_IS;
public:
  Param (int, int, std::vector<double>, std::vector<double>);
  void update_transitions(); // transition rates matrix should be altered every time
  // a new parameter is proposed or rejected, if the parameter affects the transition
  // rates. currently include lambda, mu, and interaction between serotypes
  void next_block ();
  void calc_expm(bool, int, arma::mat &, double);
  void get_rand_frailty (int num);
  void initial_calc(Data);
  double calc_llik (Data);
  double calc_lprior (int) const;
  double calc_lprior () const;
  double uni_propose (double, double, int) const;
  void alter_param (bool);
  void mcmc_move (Data, bool, double);
  double operator [] (int);
  void initialize_file(std::string);
  void print_to_file(std::string, int, arma::mat &);
};


Param::Param (int ntypes, int ntot, std::vector <double> inputs,
              std::vector <double> input_var) {
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
  proposal_sd = input_var;
  llik = 0.0;
  lprior = 0.0;
  transitions = arma::zeros(n_tot+1, n_tot+1);
  stationary_prev = arma::zeros(n_tot+1, n_tot+1);
  transitions_t = arma::zeros(n_tot+1, n_tot+1);
  update_transitions();
}
void Param::update_transitions () {
  double tot_rate = 0.0;
  for (int i=1; i<=n_tot; i++) {
    transitions(i) = mu[i-1];
    transitions(0, i) = lambda[i-1];
    tot_rate += transitions(0, i);
  }
  transitions(0) = -tot_rate;
  for (int row_i=1; row_i<=n_tot; row_i++) {
    tot_rate = transitions(row_i);
    for (int col_i=1; col_i<=n_tot; col_i++) {
      if (row_i!=col_i) {
        transitions(row_i, col_i) = lambda[col_i-1]*interaction;
        tot_rate += transitions(row_i, col_i);
      }
    }
    transitions(row_i, row_i) = -tot_rate;
  }
  stationary_prev = arma::expmat(transitions*50000); // stationary prevalence
}

void Param::next_block() {
  block_ptr++;
  if (block_ptr >= n_blocks) block_ptr = 0;
}

void Param::calc_expm(bool vacc, int ind, arma::mat& matrix_to_change, double multiplier) {
  double tot_rate = 0.0;
  for (int i=0; i<n_tot; i++) {
    matrix_to_change(0, i+1) = transitions(0, i+1)*multiplier;
    if (vacc) {
      matrix_to_change(0, i+1) *= thetaSI[i]*ind_frailty_SI[ind];
      matrix_to_change(i+1, 0) *= thetaIS[i]*ind_frailty_IS[ind];
    }
    tot_rate += matrix_to_change(0, i+1);
  }
  matrix_to_change(0, 0) = -tot_rate;
  for (int row_i=1; row_i<=n_tot; row_i++) {
    tot_rate = matrix_to_change(row_i, 0);
    for (int col_i=1; col_i<=n_tot; col_i++) {
      if (row_i!=col_i) {
        matrix_to_change(row_i, col_i) = transitions(row_i, col_i)*multiplier*interaction;
        if (vacc) {
          matrix_to_change(row_i, col_i) *= thetaSI[col_i-1]*ind_frailty_SI[ind];
        }
      }
    }
    matrix_to_change(row_i, row_i) = -tot_rate;
  }
}

double Param::calc_llik (Data data) {
  double newllik = 0.0;
  for (int i=0; i<data["n_vacc"]; i++) {
    // calculate the probability of first swab results, i.e. at stationarity
    calc_expm(true, i, stationary_prev, 50000.0);
    newllik += log(stationary_prev(0, (int)data.get_swab_v(i, 0)));
    for (int time_step=1; time_step<data["n_swabs"]; time_step++) {
      calc_expm(true, i, transitions_t, data.get_swab_time(time_step));
      newllik += log(transitions_t(data.get_swab_v(i,time_step-1),data.get_swab_v(i,time_step)));
    }
  }
  for (int i=0; i<data["n_nvacc"]; i++) {
    // calculate the probability of first swab results, i.e. at stationarity
    calc_expm(true, i, stationary_prev, 50000.0);
    newllik += log(stationary_prev(0, (int)data.get_swab_nv(i, 0)));
    for (int time_step=1; time_step<data["n_swabs"]; time_step++) {
      calc_expm(true, i, transitions_t, data.get_swab_time(time_step));
      newllik += log(transitions_t(data.get_swab_nv(i,time_step-1),data.get_swab_nv(i,time_step)));
    }
  }
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
      if (reject) {
        lambda = tempparam;
        update_transitions();
        return;
      }
      tempparam = lambda;
      for (int i=0; i<n_tot; i++) {
        lambda[i] = uni_propose(tempparam[i], proposal_sd[block_ptr], ntries);
      }
      update_transitions();
    }
    case 2: {
      if (reject) {
        mu = tempparam;
        update_transitions();
        return;
      }
      tempparam = mu;
      for (int i=0; i<n_tot; i++) {
        mu[i] = uni_propose(tempparam[i], proposal_sd[block_ptr], ntries);
      }
      update_transitions();
    }
    case 3: {
      if (reject) thetaSI = tempparam; return;
      tempparam = thetaSI;
      for (int i=0; i<n_vtypes; i++) {
        thetaSI[i] = uni_propose(tempparam[i], proposal_sd[block_ptr], ntries);
      }
    }
    case 4: {
      if (reject) thetaIS = tempparam; return;
      tempparam = thetaIS;
      for (int i=0; i<n_vtypes; i++) {
        thetaIS[i] = uni_propose(tempparam[i], proposal_sd[block_ptr], ntries);
      }
    }
    case 5: {
      if (reject) p0 = tempparam; return;
      tempparam = p0;
      for (int i=0; i<n_vtypes; i++) {
        p0[i] = uni_propose(tempparam[i], proposal_sd[block_ptr], ntries);
      }
    }
    case 6: {
      if (reject) frailtySI = tempparam[0]; return;
      tempparam[0] = frailtySI;
      frailtySI = uni_propose(tempparam[0], proposal_sd[block_ptr], ntries);
    }
    case 7: {
      if (reject) frailtyIS = tempparam[0]; return;
      tempparam[0] = frailtyIS;
      frailtyIS = uni_propose(tempparam[0], proposal_sd[block_ptr], ntries);
    }
    case 8: {
      if (reject) {
        interaction = tempparam[0];
        update_transitions();
        return;
      }
      tempparam[0] = interaction;
      interaction = uni_propose(tempparam[0], proposal_sd[block_ptr], ntries);
      update_transitions();
    }
  }
}

void Param::mcmc_move (Data data, bool adapt, double optimal_adapt) {
  // Code to propose
  alter_param (false);
  get_rand_frailty(data["n_vacc"]);
  double newllik = calc_llik(data);
  double newlprior = calc_lprior();
  double z = R::unif_rand();
  bool reject = log(z) > (newllik+newlprior-llik-lprior);
  if (reject) {
    alter_param (true);
    rejected[block_ptr]++;
  } else {
    accepted[block_ptr]++;
    llik = newllik;
    lprior = newlprior;
  }
  if (adapt) {
    double change = exp(0.999/2.0*(proposal_sd[block_ptr]-optimal_adapt));
    proposal_sd[block_ptr] *= change;
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

void Param::get_rand_frailty (int num) {
  ind_frailty_SI.clear();
  ind_frailty_IS.clear();
  for (int i=0; i<num; i++) {
    ind_frailty_SI.push_back(R::rgamma(1.0/frailtySI, frailtySI));
    ind_frailty_IS.push_back(R::rgamma(1.0/frailtyIS, frailtyIS));
  }
}

void Param::initial_calc(Data data) {
  get_rand_frailty(data["n_vacc"]);
  llik = calc_llik(data);
  lprior = calc_lprior();
}

void Param::initialize_file (std::string filename) {
  std::ofstream output_file;
  output_file.open(filename);
  output_file << "state\tposterior\tlikelihood\tprior";
  for (int i=0; i<n_tot; i++) output_file << "\tlambda" << i;
  for (int i=0; i<n_tot; i++) output_file << "\tmu" << i;
  for (int i=0; i<n_vtypes; i++) output_file << "\tthetaSI" << i;
  for (int i=0; i<n_vtypes; i++) output_file << "\tthetaIS" << i;
  for (int i=0; i<n_vtypes; i++) output_file << "\tp0" << i;
  output_file <<"\tfrailtySI\tfrailtyIS\tinteraction" << std::endl;
  output_file.close();
}

void Param::print_to_file (std::string filename, int iter, arma::mat &results_mat) {
  results_mat(iter, 0) = 0;
  results_mat(iter, 1) = llik+lprior;
  results_mat(iter, 2) = llik;
  results_mat(iter, 3) = lprior;
  std::ofstream output_file;
  output_file.open(filename, std::ofstream::app);
  output_file << iter << "\t" << llik+lprior << "\t" << llik << "\t" << lprior;
  for (int i=0; i<(n_tot*2+n_vtypes*3+3); i++) {
    output_file << "\t" << (*this)[i];
    results_mat(0, i+3) = (*this)[i];
    
  }
  output_file << std::endl;
  output_file.close();
}

bool adapt_this_iter (int iter, int adapt_every, int adapt_until, int tot_param_blocks) {
  bool adapt = (iter/tot_param_blocks)%adapt_every == 0;
  adapt = adapt & (iter*tot_param_blocks < adapt_until);
  return (adapt);
}

// [[Rcpp::export]]
arma::mat mcmc(int vaccN, int unvaccN, int n_vt, int n_nvt, int total_params,
               SEXP params, // vector of parameters
               SEXP params_sd, // vector of variances for proposal distribution
               SEXP swab_data_v, // nrow=total number of vaccinated, ncol=swabs
               SEXP swab_data_nv, // nrow=total number of un-accinated, ncol=swabs
               SEXP ab_data, // nrow=total number of vaccinated, ncol=serotypes in a vaccine
               SEXP swab_times, // timing of swabs
               SEXP mcmc_options, // options for MCMC algorithm
               SEXP filename_sexp // filename for output
                 ) {
  std::vector <double> params_vec = Rcpp::as<std::vector<double> >(params);
  std::vector <double> params_sd_vec = Rcpp::as<std::vector<double> >(params_sd);
  std::vector <double> mcmc_options_vec = Rcpp::as<std::vector<double> >(mcmc_options);
  std::string filename = Rcpp::as<std::string>(filename_sexp);
  Param parameters (n_vt, n_vt+n_nvt, params_vec, Rcpp::as<std::vector<double> >(params_sd));
  Data dataset (Rcpp::as<std::vector<double> >(swab_data_v), 
                Rcpp::as<std::vector<double> >(swab_data_nv),
                Rcpp::as<std::vector<double> >(ab_data),
                Rcpp::as<std::vector<double> >(swab_times), n_vt, n_nvt);
  std::vector <double> results (params_vec.size());
  int n_param_blocks = params_sd_vec.size();
  int niter = mcmc_options_vec[0];
  arma::mat results_mat (niter, params_vec.size()+4);
  int sample_every = mcmc_options_vec[1];
  double optimal = mcmc_options_vec[2];
  int adapt_every = mcmc_options_vec[3];
  int adapt_until = mcmc_options_vec[4];
  parameters.initialize_file(filename);
  // Calculate likelihood and prior of initial parameters
  parameters.initial_calc(dataset);
  parameters.print_to_file(filename, 0, results_mat);
  for (int iter=1; iter<niter; iter++) {
    bool adapt_bool = adapt_this_iter(iter, adapt_every, adapt_until, n_param_blocks);
    parameters.mcmc_move(dataset, adapt_bool, optimal);
    // Save parameter values
    if ((iter%sample_every)==0) {
      parameters.print_to_file(filename, 0, results_mat);
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
load("sim.data.RData")
mcmc_options <- c(niter=1000, sample_every=10, adapt_optimal=0.23, adapt_every=5, 
                  adapt_until=10)
mcmc.out <- mcmc(sim.params$N/2, sim.params$N/2, sim.params$ntypes, 
                 sim.params$ntot-sim.params$ntypes, length(sim.params.vec), 
                 unname(sim.params.vec), unname(sim.params.sd),
                 unlist(sim.data$vdata[, -1:-2]), unlist(sim.data$nvdata[, -1:-2]),
                 unlist(sim.data$abdata), sim.params$times,
                 mcmc_options, "test.txt")
*/
