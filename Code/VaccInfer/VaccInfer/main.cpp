//
//  main.cpp
//  VaccInfer
//
//  Created by Lucy Li on 2/13/17.
//  Copyright © 2017 Lucy M Li, CCDD, HPSH. All rights reserved.
//
//  compile using clang++ -L/usr/local/lib -larmadillo -O2 *.cpp > mcmc_infer
//  or in Xcode: link libarmadillo.dylib, add -larmadillo and -O2 as Other linker flags,
//  and add /usr/local/lib to library search paths

#include <iostream>
#include <stdio.h>      /* printf */
#include <math.h>       /* log */
#include <iostream>
#include <fstream>
#include "param.hpp"


bool adapt_this_iter (int iter, int adapt_every, int adapt_until, int tot_param_blocks) {
    bool adapt = (iter/tot_param_blocks)%adapt_every == 0;
    adapt = adapt & ((iter/tot_param_blocks) < adapt_until);
    if (iter < tot_param_blocks) adapt = false;
    return (adapt);
}


int main (int argc, const char * argv[]) {
    std::ifstream infile (argv[1]);
    std::string line;
    infile >> line;
    int vaccN = std::stoi(line);
    infile >> line;
    int unvaccN = std::stoi(line);
    int N = vaccN + unvaccN;
    infile >> line;
    int n_vt = std::stoi(line);
    infile >> line;
    int n_nvt = std::stoi(line);
    infile >> line;
    int total_params = std::stoi(line);
    infile >> line;
    int n_swabs = std::stoi(line);
    infile >> line;
    int niter = std::stoi(line);
    infile >> line;
    int sample_every = std::stoi(line);
    infile >> line;
    double optimal = std::stod(line);
    infile >> line;
    int adapt_every = std::stoi(line);
    infile >> line;
    int adapt_until = std::stoi(line);
    infile >> line;
    std::string filename = line; // filename for output
    infile >> line;
    int n_param_blocks = std::stoi(line);
    std::vector <int> block_sizes;
    for (int i=0; i<n_param_blocks; i++) {
        infile >> line;
        block_sizes.push_back(std::stoi(line));
    }
    std::vector <double> params_vec;
    for (int i=0; i<total_params; i++) {
        infile >> line;
        params_vec.push_back(std::stod(line));
    }
    std::vector <double> params_sd_vec;
    for (int i=0; i<n_param_blocks; i++) {
        infile >> line;
        params_sd_vec.push_back(std::stod(line));
    }
    std::vector <double> swab_data_v_vec; // nrow=total number of vaccinated, ncol=swabs
    for (int i=0; i<n_swabs*vaccN; i++) {
        infile >> line;
        swab_data_v_vec.push_back(std::stod(line));
    }
    std::vector <double> swab_data_nv_vec; // nrow=total number of un-accinated, ncol=swabs
    for (int i=0; i<n_swabs*unvaccN; i++) {
        infile >> line;
        swab_data_nv_vec.push_back(std::stod(line));
    }
    std::vector <double> ab_data_vec; // nrow=total number of vaccinated, ncol=serotypes in a vaccine
    for (int i=0; i<vaccN*n_vt; i++) {
        infile >> line;
        double value = std::stod(line);
        ab_data_vec.push_back(value);
    }
    std::vector <double> swab_times_vec; // timing of swabs
    for (int i=0; i<n_swabs; i++) {
        infile >> line;
        swab_times_vec.push_back(std::stod(line));
    }
    
    Param parameters (n_vt, n_vt+n_nvt, total_params, params_vec, params_sd_vec);
    Data dataset (swab_data_v_vec, swab_data_nv_vec, ab_data_vec, swab_times_vec, n_vt, n_nvt, vaccN, unvaccN);
    arma::mat results_mat (niter/sample_every, total_params+4);
    parameters.initialize_file(filename);
    // // Calculate likelihood and prior of initial parameters
    parameters.initial_calc(dataset);
    parameters.print_to_file(filename, 0, results_mat);
    bool adapt_bool;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    for (int iter=1; iter<niter; iter++) {
        adapt_bool = adapt_this_iter(iter, adapt_every, adapt_until, n_param_blocks);
        parameters.mcmc_move(dataset, adapt_bool, optimal);
        // Save parameter values
        if ((iter%sample_every)==0) {
            parameters.print_to_file(filename, iter/sample_every, results_mat);
            end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end-start;
            std::cout << iter << ") " << elapsed_seconds.count() << std::endl;
            start=end;
        }
    }
    return 0;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
 mcmc_options <- c(niter=5000, sample_every=100, adapt_optimal=0.23, adapt_every=10,
 adapt_until=100)
 mcmc_infer(sim.params$N/2, sim.params$N/2, sim.params$ntypes,
 sim.params$ntot-sim.params$ntypes, length(sim.params.vec),
 unname(sim.params.vec), unname(sim.params.sd),
 unlist(sim.data$vdata[, -1:-2]), unlist(sim.data$nvdata[, -1:-2]),
 c(sim.data$abdata), sim.params$times,
 mcmc_options, "test.txt")
 */

