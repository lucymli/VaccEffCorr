//
//  param.hpp
//  VaccInfer
//
//  Created by Lucy Li on 2/13/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#ifndef param_hpp
#define param_hpp

#include <iostream>
#include <cmath>
#include <armadillo>
//#include <omp.h>
#include "data.hpp"
#include "distributions.hpp"

void set_diag_as_negrowsum (arma::mat &);

class Param {
public:
    double SMALLEST_NUMBER;
    int n_vtypes, n_tot, n_params, param_index, n_ind; // counts
    std::vector <std::string> params_names;
    std::vector <double> params;
    std::vector <std::string> params_trans;
    std::vector <double> params_min, params_max;
    std::vector <std::string> params_prior;
    std::vector <double> params_prior1, params_prior2;
    std::vector <double> params_sd;
    double tempparam;
    std::vector <double> accepted, rejected;
    double llik, new_llik;
    double lprior, new_lprior;
    std::string output_file_name;
    Param ();
    Param (std::string);
    void print_params_to_screen() const;
    double transform(std::string, double, bool) const;
    void next();
    void propose();
    void accept_reject();
    double calc_lprior(bool);
    void initialize_file();
    void print_to_file(int);
    double prediction_func (double, double, double);
    void predict_lambda (arma::mat &, Data &, int, bool);
    void fill_rates (arma::mat &);
};



#endif /* param_hpp */
