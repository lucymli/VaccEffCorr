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
//    std::cout << "new_llik: " << new_llik;
//    std::cout << " new_lprior: " << new_lprior;
//    std::cout << " llik: " << llik;
//    std::cout << " lprior: " << lprior << std::endl;
    double z = log(runif());
    if (z < ratio) accept_this = true;
    if (accept_this) {
        llik = new_llik;
        lprior = new_lprior;
        accepted[block_ptr] += 1;
    }
    else {
        params.swap(tempparam);
        inferred_risk.swap(inferred_risk_temp);
        rejected[block_ptr] += 1;
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


double Param::get_param(int block_i, int i) const {
    return (params[block_starts[block_i]+i]);
}

double Param::get_inferred_risk(int ind_i, int i) const {
    double risk = inferred_risk[i*n_ind+ind_i];
    return (risk);
}

void Param::set_inferred_risk(double x, int ind_i, int i) {
    if (x > 1.0) x = 1.0;
    if (x < 0.0) x = 0.0;
    inferred_risk[i * n_ind + ind_i] = x;
}

double Param::sum_risk_by_ind(int ind_i) const {
    double risk=0.0;
    for (int i=0; i<n_vtypes; i++) {
        risk += inferred_risk[i*n_ind + ind_i];
    }
    return (risk);
}

double Param::sum_risk_by_type(int type_i) const {
    double risk=0.0;
    for (int i=0; i<n_ind; i++) {
        risk += inferred_risk[type_i*n_ind+i];
    }
    return (risk);
}

double Param::mean_risk_by_ind(int ind_i) const {
    double mean_risk = this->sum_risk_by_ind(ind_i);
    mean_risk /= (double) n_vtypes;
    return (mean_risk);
}

double Param::mean_risk_by_type(int type_i) const {
    double mean_risk = this->sum_risk_by_type(type_i);
    mean_risk /= (double) n_ind;
    return (mean_risk);
}


