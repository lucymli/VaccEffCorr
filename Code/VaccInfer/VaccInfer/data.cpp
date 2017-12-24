//
//  data.cpp
//  VaccInfer
//
//  Created by Lucy Li on 2/13/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#include <iostream>
#include "data.hpp"

Data::Data (int num_ind, int num_times, std::vector <double> time_points, int num_predictors) {
    // Generates an empty Data class
    for (int i=0; i<(num_ind*num_times); i++) carriage.push_back(0);
    for (int i=0; i<(num_ind*num_predictors); i++) {
        predictors.push_back(0);
        predictor_map.push_back(i%num_predictors);
    }
}

double Data::get_carriage(int ind_i, int time_i) const {
    return (carriage[time_i * n_ind + ind_i]);
}

void Data::set_carriage(int ind_i, int time_i, double val) {
    cariage[time_i*n_ind + ind_i] = val;
}

double Data::get_metadata(int ind_i, int i) const {
    return (metadata[i*n_ind + ind_i]);
}

double Data::get_predictor(int ind_i, int predict_i) const {
    return (predictors[predict_i * n_ind + ind_i]);
}

void Data::write_metadata_corr (int iter) const {
    std::ofstream outputfile;
    if (iter==0) {
        outputfile.open(metadata_corr_file);
        for (int i=0; i<n_covariates; i++) outputfile << "slope" << i << " intercept" << i << " R2" << i;
        outputfile << std::endl;
    }
    else {
        outputfile.open(metadata_corr_file, std::ios::app);
    }
    for (int i=0; i<metadata_corr.size(); i++) outputfile << metadata_corr[i] << " ";
    outputfile << std::endl;
    outputfile.close();
}

void Data::calc_mean_predictors() {
    for (int i=0; i<n_ind; i++) {
        mean_predictors.push_back(std::accumulate(predictors.begin()+i*n_predictors, predictors.begin()+(i+1)*n_predictors, 0.0));
        mean_predictors[n_ind] /= (double) n_predictors;
    }
}

