//
//  param.cpp
//  VaccInfer
//
//  Created by Lucy Li on 2/13/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//


#include "param.hpp"

double STATIONARY_TIME = 300.0;


void print_matrix (arma::mat & mat_to_output, std::string filename, int mat_dim) {
    std::ofstream outputfile;
    outputfile.open(filename);
    for (int row_i=0; row_i<=mat_dim; row_i++) {
        outputfile << mat_to_output(row_i, 0);
        for (int col_i=1; col_i<=mat_dim; col_i++) {
            outputfile << "\t" << mat_to_output(row_i, col_i);
        }
        outputfile << std::endl;
    }
    outputfile.close();
}


void Param::next () {
    param_index++;
    if (param_index==n_params) param_index = 0;
}

void Param::accept_reject () {
    bool accept_this = false;
    double ratio = (new_llik+new_lprior)-(llik+lprior);
    double z = log(runif());
    if (z < ratio) accept_this = true;
    if (accept_this) {
        llik = new_llik;
        lprior = new_lprior;
        accepted[param_ptr] += 1;
    }
    else {
        params[param_ptr] = tempparam;
        rejected[param_ptr] += 1;
    }
    new_llik = SMALLEST_NUMBER;
    new_lprior = SMALLEST_NUMBER;
}

void Param::print_to_file (int iter) {
    std::ofstream output_file;
    output_file.open(output_file_name, std::ofstream::app);
    output_file << iter << "\t" << llik+lprior << "\t" << llik << "\t" << lprior;
    for (int i=0; i<n_params; i++) {
        output_file << "\t" << params[i];
    }
    output_file << std::endl;
    output_file.close();
}



