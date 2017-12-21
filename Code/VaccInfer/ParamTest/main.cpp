//
//  ParamTest.cpp
//  VaccInfer
//
//  Created by Lucy Li on 12/20/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#include <stdio.h>
#include "../VaccInfer/param.hpp"


int main () {
    Param parameters("test_params.txt");
    parameters.print_params_to_screen();
    printf("====================================\n");
    printf("TEST SUCCESSFUL: Parameters loaded.\n");
    printf("====================================\n");
    double prior = parameters.calc_lprior(false);
    printf("====================================\n");
    printf("TEST SUCCESSFUL: Prior density of parameters is %f \n", prior);
    printf("====================================\n");
    parameters.param_index = parameters.n_tot;
    parameters.propose();
    std::cout << "New value proposed for " << parameters.params_names[parameters.param_index];
    std::cout << ": " << parameters.params[parameters.param_index];
    std::cout << ". Old value: " << parameters.tempparam << "." << std::endl;
    printf("====================================\n");
    double newprior = parameters.calc_lprior(false);
    printf("====================================\n");
    printf("TEST SUCCESSFUL: Prior density of newly proposed parameters is %f \n", newprior);
    std::cout << "====================================" << std::endl;
    return 0;
}

