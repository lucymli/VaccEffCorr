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


Param::Param (std::string param_file_name) {
    SMALLEST_NUMBER = -std::numeric_limits<float>::max()+100000.0;
    std::ifstream input;
    input.open(param_file_name);
    std::string name;
    param_ptr = 0;
    bool param_table_start = false;
    while (!param_table_start) {
        input >> name;
        param_table_start = name == "ParamTable";
        if (name=="output_file_name") input >> output_file_name;
        if (name=="n_vtypes") input >> n_vtypes;
        if (name=="n_tot") input >> n_tot;
        if (name=="n_params") input >> n_params;
        if (name=="n_ind") input >> n_ind;
    }
    std::string temp_str;
    double temp_double;
    for (int i=0; i<n_params; ++i) {
        input >> temp_str;
        params_names.push_back(temp_str);
        input >> temp_double;
        params.push_back(temp_double);
        params_sd.push_back(temp_double*0.5);
        input >> temp_str;
        params_trans.push_back(temp_str);
        input >> temp_double;
        params_min.push_back(temp_double);
        input >> temp_double;
        params_max.push_back(temp_double);
        input >> temp_str;
        params_prior.push_back(temp_str);
        input >> temp_double;
        params_prior1.push_back(temp_double);
        input >> temp_double;
        params_prior2.push_back(temp_double);
    }
    input.close();
    if (params.size()!=n_params) {
        std::cout << "ERROR: Expected " << n_params << " parameters but read in " << params.size() << " parameters." << std::endl;
    }
    llik = 0.0;
    new_llik = 0.0;
    lprior = 0.0;
    new_lprior = 0.0;
    param_index = 0;
    accepted.resize(n_params, 0);
    rejected.resize(n_params, 0);
}

void Param::print_params_to_screen() const {
    std::cout << "Parameter table" << std::endl;
    std::cout << "Name\tValue\tTransformation\tMin\tMax\tPrior\tPriorParam1\tPriorParam2" <<std::endl;
    for (int i=0; i<params.size(); i++) {
        std::cout << params_names[i] << "\t";
        std::cout << params[i] << "\t";
        std::cout << params_trans[i] << "\t";
        std::cout << params_min[i] << "\t";
        std::cout << params_max[i] << "\t";
        std::cout << params_prior[i] << "\t";
        std::cout << params_prior1[i] << "\t";
        std::cout << params_prior2[i] << std::endl;
    }
}

double Param::transform(std::string transformation, double value, bool reverse=false) const {
    if (transformation == "1") return (value);
    if (transformation == "inv") return (1.0/value);
    if (transformation == "log") {
        if (reverse) return (exp(value));
        return (std::log(value));
    }
    if (transformation == "log10") {
        if (reverse) return (std::pow(10, value));
        return (std::log10(value));
    }
    else return (value);
}

double Param::calc_lprior(bool print = false) {
    double dens = 0.0;
    double param_value;
    double curr_dens = 0.0;
    for (int i=0; i<params.size(); i++) {
        param_value = transform(params_trans[i], params[i]);
        curr_dens = get_density(param_value, params_prior[i], params_prior1[i], params_prior2[i], true);
        dens += curr_dens;
        if (print) std::cout << params_names[i] << ": " << curr_dens << "; cumulative: " << dens << std::endl;
    }
    if (dens < SMALLEST_NUMBER) dens = SMALLEST_NUMBER;
    return (dens);
}

void Param::propose() {
    int ntries = 50;
    tempparam = params[param_index];
    double transformed_val = transform(params_trans[param_index], tempparam);
    double rnum = rnorm(transformed_val, params_sd[param_index], params_min[param_index], params_max[param_index], ntries);
    params[param_index] = transform(params_trans[param_index], rnum, true);
}

void Param::initialize_file () {
    std::ofstream output_file;
    output_file.open(output_file_name);
    output_file << "state\tposterior\tlikelihood\tprior";
    for (int i=0; i<n_params; i++) output_file << "\t" << params_names[i];
    output_file << std::endl;
    output_file.close();
}

