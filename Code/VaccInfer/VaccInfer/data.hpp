//
//  data.hpp
//  VaccInfer
//
//  Created by Lucy Li on 2/13/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#ifndef data_hpp
#define data_hpp

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <numeric> // std::accumulate

class Data {
public:
    Data ();
    int n_ind;
    int n_time;
    int n_tot;
    int n_vtypes;
    int n_covariates;
    int n_predictors;
    std::string metadata_corr_file;
    std::vector <double> carriage; // nrow=total number of individuals, ncol=number of swabs
    std::vector <double> metadata; // nrow=total number of individuals, ncol=number of covariates
    std::vector <double> metadata_corr; //row=3, ncol=number of covariates; intercept, slope and R2
    std::vector <double> predictors; // nrow=total number of individuals, ncol=number of predictors of carriage (serotype-specific antibody levels)
    std::vector <double> mean_predictors; //nrow=total number of individuals, ncol=number of predictors of carriage
    std::vector <double> times;
    double get_carriage(int, int) const;
    double get_metadata(int, int) const;
    double get_predictor(int, int) const;
    void calc_mean_predictors();
    void write_metadata_corr (int) const;
};

#endif /* data_hpp */
