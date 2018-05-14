//
//  mcmc.cpp
//  VaccInfer
//
//  Created by Lucy Li on 5/16/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#include "mcmc.hpp"



bool adapt_this_iter (MCMC mcmc) {
    if (mcmc.iter>=mcmc.adapt_until) return (false);
    if ((mcmc.iter%mcmc.adapt_every)!=0) return (false);
    if (mcmc.iter==0) return (false);
    return (true);
}

void adapt (Param &parameters, MCMC mcmc) {
    bool adapt_bool = adapt_this_iter(mcmc);
    if (adapt_bool) {
        for (int i=0; i<parameters.n_params; i++) {
            double acceptance_rate = (double)parameters.accepted[i] / (double)(parameters.accepted[i]+parameters.rejected[i]);
            double change = exp(0.999/2.0*(acceptance_rate-mcmc.adapt_optimal));
            parameters.params_sd[i] *= change;
        }
        std::fill(parameters.accepted.begin(), parameters.accepted.end(), 0.0);
        std::fill(parameters.rejected.begin(), parameters.rejected.end(), 0.0);
    }
}


void MCMC::run_mcmc(Param parameters, Data data) {
    parameters.initialize_file();
    // // Calculate likelihood and prior of initial parameters
    parameters.lprior = parameters.calc_lprior(false);
    if (parameters.lprior > parameters.SMALLEST_NUMBER) parameters.llik = calc_llik(parameters, data, use_mean_ab);
    parameters.print_to_file(0);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    parameters.next();
    for (int iter=1; iter<niter; iter++) {
        parameters.propose();
        parameters.new_lprior = parameters.calc_lprior(false);
        if (parameters.new_lprior > parameters.SMALLEST_NUMBER) parameters.new_llik = calc_llik(parameters, data, use_mean_ab);
        parameters.accept_reject();
        adapt(parameters, *this);
        // Save parameter values
        if ((iter%sample_every)==0) {
            parameters.print_to_file(iter);
            end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end-start;
            std::cout << iter << ") " << elapsed_seconds.count() << std::endl;
            start=end;
        }
        parameters.next();
    }
}


MCMC::MCMC() {
    iter = 0;
    niter = 100000;
    sample_every = 100;
    adapt_optimal = 0.23;
    adapt_every = 10;
    adapt_until = 100;
    use_mean_ab = false;
}

MCMC::MCMC(std::string input_file) {
    std::ifstream input;
    input.open(input_file);
    std::string param_name;
    iter = 0;
    while (input) {
        input >> param_name;
        if (param_name=="niter") input >> niter;
        if (param_name=="sample_every") input >> sample_every;
        if (param_name=="adapt_optimal") input >> adapt_optimal;
        if (param_name=="adapt_every") input >> adapt_every;
        if (param_name=="adapt_until") input >> adapt_until;
        if (param_name=="use_mean_ab") input >> use_mean_ab;
    }
}

