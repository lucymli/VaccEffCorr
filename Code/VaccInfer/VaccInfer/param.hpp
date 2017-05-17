//
//  param.hpp
//  VaccInfer
//
//  Created by Lucy Li on 2/13/17.
//  Copyright © 2017 Lucy M Li, CCDD, HPSH. All rights reserved.
//

#ifndef param_hpp
#define param_hpp

#include <stdio.h>
#include <cmath>
#include <armadillo>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions.hpp>
#include <omp.h>
#include "data.hpp"

class Param {
public:
    Param ();
    std::vector <int> num_per_block;
    int n_vtypes, n_nvtypes, n_tot, n_params, n_blocks, block_ptr, param_ptr; // counts
    std::vector <double> tempparam, params;
    std::vector <double> accepted;
    double llik, new_llik;
    double lprior, new_lprior;
    
    arma::mat transitions, stationary_prev, transitions_t;
    std::string output_file_name;
    void update_transitions(); // transition rates matrix should be altered every time
    // a new parameter is proposed or rejected, if the parameter affects the transition
    // rates. currently include lambda, mu, and interaction between serotypes
    void next_block ();
//    void calc_expm(bool, int, Data, arma::mat &, arma::mat &, double);
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

#endif /* param_hpp */
