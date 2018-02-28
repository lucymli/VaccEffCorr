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
#include "param.hpp"
#include "data.hpp"



double calc_llik (Param, Data, bool);

void get_multinom_prob (Param, Data, double, arma::mat &, bool);


#endif /* likelihood_hpp */
