#include "Param.h"

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
  transitions
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
