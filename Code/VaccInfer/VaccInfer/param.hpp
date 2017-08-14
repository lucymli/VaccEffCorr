//
//  param.hpp
//  VaccInfer
//
//  Created by Lucy Li on 2/13/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#ifndef param_hpp
#define param_hpp

#include <stdio.h>
#include <cmath>
#include <armadillo>
#include "data.hpp"
#include "distributions.hpp"

class Param {
public:
    double SMALLEST_NUMBER;
    std::vector <int> num_per_block;
    std::vector <int> block_starts;
    int n_vtypes, n_tot, n_params, n_blocks, block_ptr, param_ptr, n_ind; // counts
    std::vector <double> tempparam, params, params_sd;
    std::vector <double> accepted, rejected;
    double llik, new_llik;
    double lprior, new_lprior;
    std::string output_file_name;
    std::string inferred_risk_file;
    std::vector <double> inferred_risk; // nrow=total number of individuals, ncol=number of serotypes
    std::vector <double> inferred_risk_temp;
    Param ();
    void next();
    void dprior();
    void rprior();
    void propose();
    void accept_reject();
    double calc_lprior();
    void initialize_file();
    void print_to_file(int);
    double get_param(int, int) const;
    double get_inferred_risk(int, int) const;
    void set_inferred_risk(double, int, int);
    double sum_risk_by_ind(int) const;
    void write_inferred_risk (int) const;
    double sum_risk_by_type(int) const;
    double mean_risk_by_ind(int) const;
    double mean_risk_by_type(int) const;
};

#endif /* param_hpp */
