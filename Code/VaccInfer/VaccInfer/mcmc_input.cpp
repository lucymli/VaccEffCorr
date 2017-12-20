//
//  mcmc.cpp
//  VaccInfer
//
//  Created by Lucy Li on 5/16/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//


#include "mcmc.hpp"


MCMC::MCMC() {
    iter = 0;
    niter = 100000;
    sample_every = 100;
    adapt_optimal = 0.23;
    adapt_every = 10;
    adapt_until = 100;
    tot_blocks = 5;
    use_mean_ab = false;
}

MCMC::MCMC(std::string input_file) {
    std::ifstream input;
    input.open(input_file);
    std::string param_name;
    while (input) {
        input >> param_name;
        if (param_name=="iter") input >> iter;
        if (param_name=="niter") input >> niter;
        if (param_name=="sample_every") input >> sample_every;
        if (param_name=="adapt_optimal") input >> adapt_optimal;
        if (param_name=="adapt_every") input >> adapt_every;
        if (param_name=="adapt_until") input >> adapt_until;
        if (param_name=="tot_blocks") input >> tot_blocks;
        if (param_name=="use_mean_ab") input >> use_mean_ab;
    }
}
