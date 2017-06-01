//
//  data.cpp
//  VaccInfer
//
//  Created by Lucy Li on 2/13/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#include <iostream>
#include "data.hpp"


double Data::get_carriage(int ind_i, int time_i) const {
    return (carriage[time_i * n_ind + ind_i]);
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

