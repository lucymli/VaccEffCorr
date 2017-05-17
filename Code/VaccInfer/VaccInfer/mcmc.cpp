//
//  mcmc.cpp
//  VaccInfer
//
//  Created by Lucy Li on 5/16/17.
//  Copyright © 2017 Lucy M Li, CCDD, HPSH. All rights reserved.
//

#include "mcmc.hpp"


bool adapt_this_iter (MCMC mcmc) {
    if ((mcmc.iter/mcmc.tot_blocks) < mcmc.adapt_every) return (false);
    bool adapt = (mcmc.iter/mcmc.tot_blocks)%mcmc.adapt_every == 0;
    adapt = adapt & (mcmc.iter < (mcmc.adapt_until*mcmc.tot_blocks));
    return (adapt);
}

void print_parameters (Param parameters) {
    
}

double calc_lprior(Param parameters, bool log) {
    double lprior = 0.0;
    return (lprior);
}

void propose (Param &parameters) {
}

void accept_reject (Param &parameters) {
    
}


void adapt(Param &parameters, MCMC mcmc) {
    bool adapt_bool = adapt_this_iter(mcmc);
    if (adapt_bool);
}

double lm_intercept(Param &parameters, Data data) {
}


void MCMC::run_mcmc(Param parameters, Data data) {
    print_parameters(parameters);
    // // Calculate likelihood and prior of initial parameters
    parameters.llik = calc_lik(parameters, data, true);
    parameters.lprior = calc_lprior(parameters, true);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    for (int iter=1; iter<niter; iter++) {
        propose(parameters);
        parameters.new_llik = calc_lik(parameters, data, true);
        parameters.new_lprior = calc_lprior(parameters, true);
        accept_reject(parameters);
        adapt(parameters, *this);
        // Save parameter values
        if ((iter%sample_every)==0) {
            print_parameters(parameters);
            end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end-start;
            std::cout << iter << ") " << elapsed_seconds.count() << std::endl;
            start=end;
        }
    }
}
