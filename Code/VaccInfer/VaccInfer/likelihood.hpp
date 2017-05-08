//
//  likelihood.hpp
//  VaccInfer
//
//  Created by Lucy Li on 5/8/17.
//  Copyright © 2017 Lucy M Li, CCDD, HPSH. All rights reserved.
//

#ifndef likelihood_hpp
#define likelihood_hpp

#include <stdio.h>
#include "param.hpp"
#include "data.hpp"



double calc_lik (Param, Data, bool);

double predict_lambda (double, Param);


#endif /* likelihood_hpp */
