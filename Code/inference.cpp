#include <RcppArmadillo.h>
#include <stdio.h>      /* printf */
#include <math.h>       /* log */
#include <iostream>
#include <fstream>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double SMALLEST_NUMBER = -std::numeric_limits<float>::max()+100000.0;
double STATIONARY_TIME = 365000.0;

void print_matrix (arma::mat & mat_to_output, std::string filename, int mat_dim) {
  std::ofstream outputfile;
  outputfile.open(filename);
  for (int row_i=0; row_i<=mat_dim; row_i++) {
    outputfile << mat_to_output(row_i, 0);
    for (int col_i=1; col_i<=mat_dim; col_i++) {
      outputfile << "\t" << mat_to_output(row_i, col_i);
    }
    outputfile << std::endl;
  }
  outputfile.close();
}

class Data {
  std::vector <double> swab_data_v; // nrow=total number of vaccinated, ncol=swabs
  std::vector <double> swab_data_nv; // nrow=total number of vaccinated, ncol=swabs
  std::vector <double> swab_times; // timing of swabs
  std::vector <double> ab_data; // nrow=total number of vaccinated, ncol=serotypes in a vaccine
  std::vector <double> lambda_reduction; // nrow=total number of vaccinated, ncol=serotypes in a vaccine
  int n_vacc, n_nvacc, n_swabs, n_vtypes, n_nvtypes,n_tot;
  std::vector<double> mean_ab, max_ab;
public:
  Data (std::vector<double>, std::vector<double>, std::vector<double>,
        std::vector <double>, int, int, int, int);
  int operator[] (std::string) const;
  double get_swab_time (int);
  double get_swab_v (int, int);
  double get_swab_nv (int, int);
  double get_ab (int, int);
  double get_mean_ab (int);
  double get_max_ab (int);
};


Data::Data (std::vector<double> vdata, std::vector<double>nvdata, 
            std::vector<double> abdata, std::vector <double>timings,
            int vtypes, int nvtypes, int nvacc, int nnvacc) {
  for (auto i=vdata.begin(); i<vdata.end(); i++) swab_data_v.push_back(*i);
  for (auto i=nvdata.begin(); i<nvdata.end(); i++) swab_data_nv.push_back(*i);
  for (auto i=abdata.begin(); i<abdata.end(); i++) ab_data.push_back(*i);
  swab_times.push_back(timings[0]);
  for (int i=1; i<timings.size(); i++) swab_times.push_back(timings[i]-timings[i-1]);
  n_swabs = timings.size();
  n_vacc = nvacc;
  n_nvacc = nnvacc;
  n_vtypes = vtypes;
  n_nvtypes = nvtypes;
  n_tot = n_vtypes + n_nvtypes;
  for (int i=0; i<n_vtypes; i++) {
    auto max_ele = max_element(abdata.begin()+i*n_vacc, abdata.begin()+(i+1)*n_vacc);
    max_ab.push_back(*max_ele);
    double mean_val = 0.0;
    for (int j=0; j<n_vacc; j++) {
      mean_val += max_ab.back()-abdata[i*n_vacc+j];
      abdata[i*n_vacc+j] = max_ab.back()-abdata[i*n_vacc+j];
    }
    mean_ab.push_back(mean_val/(double)n_vacc);
  }
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

double Data::get_ab (int ind_i, int type_i) {
  return(ab_data[type_i*n_vacc+ind_i]);
}

double Data::get_mean_ab (int type_i) {
  return(mean_ab[type_i]);
}

double Data::get_max_ab (int type_i) {
  return(max_ab[type_i]);
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
  arma::mat transitions, stationary_prev, transitions_t;
  std::vector <double> ind_frailty_SI, ind_frailty_IS;
public:
  Param (int, int, std::vector<double>, std::vector<double>);
  void update_transitions(); // transition rates matrix should be altered every time
  // a new parameter is proposed or rejected, if the parameter affects the transition
  // rates. currently include lambda, mu, and interaction between serotypes
  void next_block ();
  void calc_expm(bool, int, Data, arma::mat &, double);
  void get_rand_frailty (Data &);
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
  /* Check the parameters have been read in correctly */
  // for (auto i=lambda.begin(); i<lambda.end(); i++) Rcout << "lambda " << *i << std::endl;
  // for (auto i=mu.begin(); i<mu.end(); i++) Rcout << "mu " << *i << std::endl;
  // for (auto i=p0.begin(); i<p0.end(); i++) Rcout << "p0 " << *i << std::endl;
  // for (auto i=thetaSI.begin(); i<thetaSI.end(); i++) Rcout << "thetaSI " << *i << std::endl;
  // for (auto i=thetaIS.begin(); i<thetaIS.end(); i++) Rcout << "thetaIS " << *i << std::endl;
  // Rcout << "frailtySI: " << frailtySI << " frailtyIS: " << frailtyIS << " interaction: " << interaction << std::endl;
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
    transitions(i, 0) = mu[i-1];
    transitions(0, i) = lambda[i-1];
    tot_rate += transitions(0, i);
  }
  transitions(0, 0) = -tot_rate;
  for (int row_i=1; row_i<=n_tot; row_i++) {
    tot_rate = transitions(row_i, 0);
    for (int col_i=1; col_i<=n_tot; col_i++) {
      if (row_i!=col_i) {
        transitions(row_i, col_i) = lambda[col_i-1]*interaction;
        tot_rate += transitions(row_i, col_i);
      }
    }
    transitions(row_i, row_i) = -tot_rate;
  }
}

void Param::next_block() {
  block_ptr++;
  if (block_ptr >= n_blocks) block_ptr = 0;
}

void Param::calc_expm(bool vacc, int ind, Data data, arma::mat& original_matrix, double multiplier) {
  int n_vacc = data["n_vacc"];
  arma::mat matrix_to_change(original_matrix);
  double tot_rate = 0.0;
  Rcout << ind << "__" << std::endl;
  for (int i=0; i<n_tot; i++) {
    double full_immunity = (double) (ind <= (n_vacc*p0[i]));
    matrix_to_change(0, i+1) = transitions(0, i+1)*multiplier;
    matrix_to_change(i+1, 0) = transitions(i+1, 0)*multiplier;
    if (vacc & (i<=n_vtypes)) {
      matrix_to_change(0, i+1) *= thetaSI[i]*full_immunity*(data.get_max_ab(i-1)-data.get_ab(ind, i-1))/data.get_max_ab(i-1);//*ind_frailty_SI[(i-1)*n_vacc+ind]
      matrix_to_change(i+1, 0) *= thetaIS[i]*full_immunity;//*ind_frailty_IS[ind];
    }
    if (matrix_to_change(0, i+1)<0.0) matrix_to_change(0, i+1) = 0.0;
    tot_rate += matrix_to_change(0, i+1);
  }
  matrix_to_change(0, 0) = -tot_rate;
  for (int row_i=1; row_i<=n_tot; row_i++) {
    tot_rate = matrix_to_change(row_i, 0);
    double full_immunity=1.0;
    if (row_i<=n_vtypes) full_immunity = (double) (ind <= (n_vacc*p0[row_i-1]));
    for (int col_i=1; col_i<=n_tot; col_i++) {
      if (row_i!=col_i) {
        matrix_to_change(row_i, col_i) = transitions(row_i, col_i)*multiplier*interaction;
        if (vacc & (col_i<=n_vtypes)) {
          matrix_to_change(row_i, col_i) *= thetaSI[col_i-1]*full_immunity*(data.get_max_ab(col_i-1)-data.get_ab(ind, col_i-1))/data.get_max_ab(col_i-1);//*ind_frailty_SI[(col_i-1)*n_vacc+ind]
        }
        if (matrix_to_change(row_i, col_i)<0.0) matrix_to_change(row_i, col_i) = 0.0;
        tot_rate += matrix_to_change(row_i, col_i);
      }
      else matrix_to_change(row_i, col_i) = 0.0;
    }
    matrix_to_change(row_i, row_i) = -tot_rate;
  }
  try {
    original_matrix = arma::expmat(matrix_to_change);
  }
  catch(...) {
    print_matrix(matrix_to_change, "matrix.txt", n_tot);
  }
}

double Param::calc_llik (Data data) {
  double newllik = 0.0;
  int position1, position2;
  for (int i=0; i<data["n_vacc"]; i++) {
    /*For each individual, calculate the likelihood of observations*/
    position1 = (int) data.get_swab_v(i, 0);
    // calculate the probability of first swab results, i.e. at stationarity
    try{
      calc_expm(true, i, data, stationary_prev, STATIONARY_TIME);
    }
    catch (...) {
      return (SMALLEST_NUMBER);
    }
    newllik += log(stationary_prev(0, position1));
    for (int time_step=1; time_step<data["n_swabs"]; time_step++) {
      try{
        calc_expm(true, i, data, transitions_t, data.get_swab_time(time_step));
      }
      catch (...) {
        return (SMALLEST_NUMBER);
      }
      position2 = (int) data.get_swab_v(i, time_step);
      newllik += log(transitions_t(position1,position2));
      // Rcout << "i: "<< i << " time_step: " << time_step << " lik: " << newllik << std::endl;
      position1 = position2;
    }
  }
  try{
    calc_expm(false, 0, data, stationary_prev, STATIONARY_TIME);
  }
  catch (...) {
    return (SMALLEST_NUMBER);
  }
  for (int i=0; i<data["n_nvacc"]; i++) {
    position1 = (int) data.get_swab_nv(i, 0);
    newllik += log(stationary_prev(0, position1));
    for (int time_step=1; time_step<data["n_swabs"]; time_step++) {
      try{
        calc_expm(false, i, data, transitions_t, data.get_swab_time(time_step));
      }
      catch (...) {
        return (SMALLEST_NUMBER);
      }
      position2 = (int) data.get_swab_nv(i, time_step);
      newllik += log(transitions_t(position1, position2));
      position1 = position2;
    }
  }
  if (!std::isfinite(newllik)) newllik = SMALLEST_NUMBER;
  return (newllik);
}

double Param::calc_lprior (int block_i) const {
  double newlprior = 0.0;
  if (block_i==0) for (int i=0; i<lambda.size(); i++) newlprior += R::dlnorm(lambda[i], -5.9, 0.2, 1);
  else if (block_i==1) for (auto i=mu.begin(); i<mu.end(); i++) newlprior += R::dlnorm(*i, -5.9, 0.2, 1);
  else if (block_i==2) for (auto i=thetaSI.begin(); i<thetaSI.end(); i++) newlprior += R::dunif(*i, 0, 1, 1);
  else if (block_i==3) for (auto i=thetaIS.begin(); i<thetaIS.end(); i++) newlprior += R::dunif(*i, 0.5, 2.0, 1);
  else if (block_i==4) for (auto i=p0.begin(); i<p0.end(); i++) newlprior += R::dunif(*i, 0.0, 1.0, 1);
  else if (block_i==5) newlprior += R::dexp(frailtySI, 1, 1);
  else if (block_i==6) newlprior += R::dexp(frailtyIS, 1, 1);
  else if (block_i==7) newlprior += R::dbeta(interaction, 200, 300, 1); // mean = 0.4, sd = 0.2
  if (!std::isfinite(newlprior)) newlprior = SMALLEST_NUMBER;
  return (newlprior);
}

double Param::calc_lprior () const {
  double newlprior = 0.0;
  for (int i=0; i<n_blocks; i++) {
    newlprior += calc_lprior(i);
    if (newlprior <= SMALLEST_NUMBER) return (SMALLEST_NUMBER);
  }
  return (newlprior);
}

double Param::uni_propose (double oldval, double sd, int ntries) const {
  double newval = R::rnorm(oldval, sd);
  // while (newval < 0) {
  newval = R::rnorm(oldval, sd);
  // ntries--;
  // if (ntries < 0) break;
  // }
  return (newval);
}

void Param::alter_param (bool reject) {
  int ntries = 10;
  if (block_ptr==0) {
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
  else if (block_ptr==1) {
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
  else if (block_ptr==2) {
    if (reject) thetaSI = tempparam; return;
    tempparam = thetaSI;
    for (int i=0; i<n_vtypes; i++) {
      thetaSI[i] = uni_propose(tempparam[i], proposal_sd[block_ptr], ntries);
    }
  }
  else if (block_ptr==3) {
    if (reject) thetaIS = tempparam; return;
    tempparam = thetaIS;
    for (int i=0; i<n_vtypes; i++) {
      thetaIS[i] = uni_propose(tempparam[i], proposal_sd[block_ptr], ntries);
    }
  }
  else if (block_ptr==4) {
    if (reject) p0 = tempparam; return;
    tempparam = p0;
    for (int i=0; i<n_vtypes; i++) {
      p0[i] = uni_propose(tempparam[i], proposal_sd[block_ptr], ntries);
    }
  }
  else if (block_ptr==5) {
    if (reject) frailtySI = tempparam[0]; return;
    tempparam[0] = frailtySI;
    frailtySI = uni_propose(tempparam[0], proposal_sd[block_ptr], ntries);
  }
  else if (block_ptr==6) {
    if (reject) frailtyIS = tempparam[0]; return;
    tempparam[0] = frailtyIS;
    frailtyIS = uni_propose(tempparam[0], proposal_sd[block_ptr], ntries);
  }
  else if (block_ptr==7) {
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

void Param::mcmc_move (Data data, bool adapt, double optimal_adapt) {
  if ((block_ptr == 3) | (block_ptr == 6) | (block_ptr == 7)) next_block();
  // Code to propose
  alter_param (false);
  get_rand_frailty(data);
  double newlprior = calc_lprior();
  double newllik = SMALLEST_NUMBER;
  if (newlprior > SMALLEST_NUMBER) newllik = calc_llik(data);
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
    double acceptance_rate = (double)accepted[block_ptr] / (double)(accepted[block_ptr]+rejected[block_ptr]);
    double change = exp(0.999/2.0*(acceptance_rate-optimal_adapt));
    proposal_sd[block_ptr] *= change;
    std::fill(accepted.begin(), accepted.end(), 0.0);
    std::fill(rejected.begin(), rejected.end(), 0.0);
  }
  next_block();
}


double Param::operator[](int i) {
  double value;
  if (i < n_tot) value = lambda[i];
  else if (i < n_tot*2) value = mu[i-n_tot];
  else if (i < n_tot*2+n_vtypes) value = thetaSI[i-n_tot*2];
  else if (i < n_tot*2+n_vtypes*2) value = thetaIS[i-n_tot*2-n_vtypes];
  else if (i < n_tot*2+n_vtypes*3) value = p0[i-n_tot*2-n_vtypes*2];
  else if (i == n_tot*2+n_vtypes*3) value = frailtySI;
  else if (i == n_tot*2+n_vtypes*3+1) value = frailtyIS;
  else value = interaction;
  return (value);
}

void Param::get_rand_frailty (Data & data) {
  ind_frailty_SI.clear();
  // ind_frailty_IS.clear();
  for (int type_i=0; type_i<data["n_vtypes"]; type_i++) {
    double m1 = data.get_mean_ab(type_i);
    for (int ind_i=0; ind_i<data["n_vacc"]; ind_i++) {
      // ind_frailty_SI.push_back((data.get_ab(ind_i, type_i)*sqrt(frailtySI)/sd1 + (1.0-m1)));
      ind_frailty_SI.push_back((data.get_max_ab(type_i)-data.get_ab(ind_i, type_i))/m1);
    }
    // ind_frailty_SI.push_back(R::rgamma(1.0/frailtySI, frailtySI));
    // ind_frailty_IS.push_back(R::rgamma(1.0/frailtyIS, frailtyIS));
  }
}

void Param::initial_calc(Data data) {
  // get_rand_frailty(data);
  lprior = calc_lprior();
  if (lprior > SMALLEST_NUMBER) llik = calc_llik(data);
  else llik = SMALLEST_NUMBER;
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

void Param::print_to_file (std::string filename, int iter, arma::mat &output_mat) {
  output_mat(iter, 0) = iter;
  output_mat(iter, 1) = llik+lprior;
  output_mat(iter, 2) = llik;
  output_mat(iter, 3) = lprior;
  std::ofstream output_file;
  output_file.open(filename, std::ofstream::app);
  output_file << iter << "\t" << llik+lprior << "\t" << llik << "\t" << lprior;
  for (int i=0; i<(n_tot*2+n_vtypes*3+3); i++) {
    output_file << "\t" << (*this)[i];
    output_mat(iter, i+4) = (*this)[i];
  }
  output_file << std::endl;
  output_file.close();
}

bool adapt_this_iter (int iter, int adapt_every, int adapt_until, int tot_param_blocks) {
  bool adapt = (iter/tot_param_blocks)%adapt_every == 0;
  adapt = adapt & ((iter/tot_param_blocks) < adapt_until);
  if (iter < tot_param_blocks) adapt = false;
  return (adapt);
}

// [[Rcpp::export]]
arma::mat mcmc_infer(int vaccN, int unvaccN, int n_vt, int n_nvt, int total_params,
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
  std::vector <double> swab_data_v_vec = Rcpp::as<std::vector<double> >(swab_data_v);
  std::vector <double> swab_data_nv_vec = Rcpp::as<std::vector<double> >(swab_data_nv);
  std::vector <double> ab_data_vec = Rcpp::as<std::vector<double> >(ab_data);
  std::vector <double> swab_times_vec = Rcpp::as<std::vector<double> >(swab_times);
  Data dataset (swab_data_v_vec, swab_data_nv_vec, ab_data_vec, swab_times_vec, n_vt, n_nvt, vaccN, unvaccN);
  std::vector <double> results (params_vec.size());
  int n_param_blocks = params_sd_vec.size();
  int niter = mcmc_options_vec[0];
  int sample_every = mcmc_options_vec[1];
  double optimal = mcmc_options_vec[2];
  int adapt_every = mcmc_options_vec[3];
  int adapt_until = mcmc_options_vec[4];
  arma::mat results_mat (niter/sample_every, params_vec.size()+4);
  parameters.initialize_file(filename);
  // // Calculate likelihood and prior of initial parameters
  parameters.initial_calc(dataset);
  parameters.print_to_file(filename, 0, results_mat);
  bool adapt_bool;
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  for (int iter=1; iter<niter; iter++) {
    adapt_bool = adapt_this_iter(iter, adapt_every, adapt_until, n_param_blocks);
    parameters.mcmc_move(dataset, adapt_bool, optimal);
    // Save parameter values
    if ((iter%sample_every)==0) {
      parameters.print_to_file(filename, iter/sample_every, results_mat);
      end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end-start;
      Rcout << iter << ") " << elapsed_seconds.count() << std::endl;
      start=end;
    }
  }
  return results_mat;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
setwd("~/Dropbox (VERG)/Research/Projects/Vaccine/VaccEffCorr/Code")
source("parse_data.R")
sim.data <- new.env()
load("sim.data.RData", sim.data)
plot.prev.timeseries(data.frame(apply(sim.data$sim.data$vdata[, -1:-2], 2, findInterval, c(1, sim.data$sim.params$ntypes))))
load("sim.data.RData")
mcmc_options <- c(niter=5000, sample_every=100, adapt_optimal=0.23, adapt_every=10,
                  adapt_until=100)
mcmc.out <- data.frame(mcmc_infer(sim.params$N/2, sim.params$N/2, sim.params$ntypes,
                                  sim.params$ntot-sim.params$ntypes, length(sim.params.vec),
                                  unname(sim.params.vec), unname(sim.params.sd),
                                  unlist(sim.data$vdata[, -1:-2]), unlist(sim.data$nvdata[, -1:-2]),
                                  c(sim.data$abdata), sim.params$times,
                                  mcmc_options, "test.txt"))
library(ggplot2)
selected.var <- c(1:4, 5:(4+sim.params$ntot*2+sim.params$ntypes), (5+sim.params$ntot*2+sim.params$ntypes*2):(5+sim.params$ntot*2+sim.params$ntypes*3+1), ncol(mcmc.out))
# ggplot(reshape2::melt(mcmc.out[, selected.var], id.vars="X1"), aes(x=X1, y=value)) + theme_bw() +
  # geom_line() + facet_grid(variable~., scales="free_y")
# cor(mcmc.out[, c(10:14, ncol(mcmc.out))])
*/
