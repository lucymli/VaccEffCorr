//
//  mcmc.cpp
//  VaccInfer
//
//  Created by Lucy Li on 5/16/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#include "mcmc.hpp"



bool adapt_this_iter (MCMC mcmc) {
    if ((mcmc.iter/mcmc.tot_blocks) < mcmc.adapt_every) return (false);
    bool adapt = (mcmc.iter/mcmc.tot_blocks)%mcmc.adapt_every == 0;
    adapt = adapt & (mcmc.iter < (mcmc.adapt_until*mcmc.tot_blocks));
    return (adapt);
}

void adapt (Param &parameters, MCMC mcmc) {
    bool adapt_bool = adapt_this_iter(mcmc);
    if (adapt_bool) {
        double acceptance_rate;
        acceptance_rate = (double)parameters.accepted[parameters.block_ptr] /
            (double)(parameters.accepted[parameters.block_ptr]+parameters.rejected[parameters.block_ptr]);
        double change = exp(0.999/2.0*(acceptance_rate-mcmc.adapt_optimal));
        parameters.params_sd[parameters.block_ptr] *= change;
        std::fill(parameters.accepted.begin(), parameters.accepted.end(), 0.0);
        std::fill(parameters.rejected.begin(), parameters.rejected.end(), 0.0);
    }
}

void linear_regression (Param &parameters, Data data) {
    double n = (double) data.n_ind;
    std::vector <double> y;
    for (int ind_i=0; ind_i<data.n_ind; ind_i++) y.push_back(parameters.mean_risk_by_ind(ind_i));
    double sum_y = std::accumulate(y.begin(), y.end(), 0.0);
    double mean_y = sum_y / (double) y.size();
    double intercept, slope, numerator, denominator, R2;
    double sum_x, sum_xy, sum_x2, x;
    for (int metadata_i=0; metadata_i<data.n_covariates; metadata_i++) {
        sum_x = 0.0, sum_x2 = 0.0;
        sum_xy = 0.0;
        for (int ind_i=0; ind_i<data.n_ind; ind_i++) {
            x = data.metadata[metadata_i*data.n_ind+ind_i];
            sum_x += x;
            sum_x2 += x * x;
            sum_xy += x * y[ind_i];
        }
        intercept = (sum_y*sum_x2-sum_x*sum_xy)/(n*sum_x2-sum_x*sum_x);
        slope = (n*sum_xy-sum_x*sum_y) / (n*sum_x2-sum_x*sum_x);
        numerator=0.0, denominator=0.0;
        for (int ind_i=0; ind_i<data.n_ind; ind_i++) {
            double x = data.metadata[metadata_i*data.n_ind+ind_i];
            numerator += pow(intercept + slope * x - mean_y, 2.0);
            denominator += pow(y[ind_i] - mean_y, 2.0);
        }
        R2 = 1.0 - numerator / denominator;
        data.metadata_corr[metadata_i*3+0] = intercept;
        data.metadata_corr[metadata_i*3+1] = slope;
        data.metadata_corr[metadata_i*3+2] = R2;
    }
}


void MCMC::run_mcmc(Param parameters, Data data) {
    parameters.initialize_file();
    // // Calculate likelihood and prior of initial parameters
    parameters.lprior = parameters.calc_lprior();
    if (parameters.lprior > parameters.SMALLEST_NUMBER) parameters.llik = calc_llik(parameters, data, use_mean_ab);
    parameters.print_to_file(0);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
//    parameters.write_metadata_corr(0);
    //parameters.write_inferred_risk(0);
    parameters.next();
    for (int iter=1; iter<niter; iter++) {
        parameters.propose();
        parameters.new_lprior = parameters.calc_lprior();
        if (parameters.new_lprior > parameters.SMALLEST_NUMBER) parameters.new_llik = calc_llik(parameters, data, use_mean_ab);
        parameters.accept_reject();
        adapt(parameters, *this);
        // Save parameter values
        if ((iter%sample_every)==0) {
            parameters.print_to_file(iter);
//            parameters.write_metadata_corr(iter);
            //parameters.write_inferred_risk(iter);
            end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end-start;
            std::cout << iter << ") " << elapsed_seconds.count() << std::endl;
            start=end;
        }
        parameters.next();
    }
}
