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
    std::vector <double> metadata_corr; //row=3, ncol=number of covariates; intercept, slope and R2
    std::vector <double> predictors; // nrow=total number of individuals, ncol=number of predictors of carriage (serotype-specific antibody levels)
    std::vector <double> inferred_risk; // nrow=total number of individuals, ncol=number of serotypes
    double get_carriage(int, int) const;
    double get_metadata(int, int) const;
    double get_predictor(int, int) const;
    double get_inferred_risk(int, int) const;
    double sum_risk_by_ind(int) const;
    double sum_risk_by_type(int) const;
    double mean_risk_by_ind(int) const;
    double mean_risk_by_type(int) const;
    int n_ind;
    int n_time;
    int n_types;
    int n_vacc_types;
    int n_covariates;
    int n_predictors;
};

#endif /* data_hpp */
