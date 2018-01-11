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
    Data (std::string);
    Data (int, int, int, std::vector <double>, int); // generates an empty Data class
    int n_ind;
    int n_time;
    int n_tot;
    int n_vtypes;
    int n_covariates;
    int n_predictors;
    std::string metadata_corr_file;
    std::vector <int> carriage; // nrow=total number of individuals, ncol=number of swabs
    std::vector <double> metadata; // nrow=total number of individuals, ncol=number of covariates
    std::vector <double> metadata_corr; //row=3, ncol=number of covariates; intercept, slope and R2
    std::vector <double> predictors; // nrow=total number of individuals, ncol=number of predictors of carriage (serotype-specific antibody levels)
    std::vector <int> predictor_map;
    std::vector <double> mean_predictors; //nrow=total number of individuals, ncol=number of predictors of carriage
    std::vector <double> times;
    int get_carriage(int, int) const;
    void set_carriage(int, int, int);
    double get_metadata(int, int) const;
    double get_predictor(int, int) const;
    void set_predictor(int, int, double);
    int get_predictor_index(int, int) const;
    void set_predictor_index(int ind_i, int predictor_i, int mapping);
    void calc_mean_predictors();
    void print_data_to_file(std::string);
};


#endif /* data_hpp */
