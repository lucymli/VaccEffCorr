//
//  data.hpp
//  VaccInfer
//
//  Created by Lucy Li on 2/13/17.
//  Copyright © 2017 Lucy M Li, CCDD, HPSH. All rights reserved.
//

#ifndef data_hpp
#define data_hpp

#include <stdio.h>
#include <vector>
#include <string>

class Data {
public:
    Data ();
    std::vector <double> carriage; // nrow=total number of individuals, ncol=number of swabs
    std::vector <double> metadata; // nrow=total number of individuals, ncol=number of covariates
    std::vector <double> predictors; // nrow=total number of individuals, ncol=number of predictors of carraige (serotype-specific antibody levels)
    int n_ind;
    int n_time;
    int n_types;
    int n_vacc_types;
    int n_covariates;
    int n_predictors;
};

#endif /* data_hpp */
