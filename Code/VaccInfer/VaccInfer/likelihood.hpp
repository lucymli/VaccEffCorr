//
//  likelihood.hpp
//  VaccInfer
//
//  Created by Lucy Li on 5/8/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#ifndef likelihood_hpp
#define likelihood_hpp

#include <stdio.h>
#include <armadillo>
#include <omp.h>
#include "param.hpp"
#include "data.hpp"

double calc_llik (Param &, Data&, bool);

void fill_rates (Param, arma::mat &);

arma::mat fill_rates (arma::mat, Data, Param, int);

double predict_lambda (Param&, Data, bool);


#endif /* likelihood_hpp */
