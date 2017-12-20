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
    printf("TEST SUCCESSFUL: Parameters loaded.\n");
    parameters.print_params_to_screen();
    double prior = parameters.calc_lprior(true);
    printf("TEST SUCCESSFUL: Prior density is %f \n", prior);
    return 0;
}

