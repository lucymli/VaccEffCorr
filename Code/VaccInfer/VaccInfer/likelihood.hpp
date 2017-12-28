//
//  likelihood.hpp
//  VaccInfer
//
//  Created by Lucy Li on 5/8/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#ifndef likelihood_hpp
#define likelihood_hpp

#include <iostream>
#include <armadillo>
//#include <omp.h>
#include "param.hpp"
#include "data.hpp"

void set_diag_as_negrowsum (arma::mat &);

double prediction_func (double, double, double);

void predict_lambda (arma::mat &, Param &, Data &, int, bool);

void fill_rates (Param &, arma::mat &);

double calc_llik (Param &, Data&, bool);


#endif /* likelihood_hpp */
