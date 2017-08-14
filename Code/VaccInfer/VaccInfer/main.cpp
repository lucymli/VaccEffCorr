//
//  main.cpp
//  VaccInfer
//
//  Created by Lucy Li on 2/13/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//
//  compile using clang++ -L/usr/local/lib -larmadillo -O2 *.cpp > mcmc_infer
//  or in Xcode: link libarmadillo.dylib, add -larmadillo and -O2 as Other linker flags,
//  and add /usr/local/lib to library search paths

#include <iostream>
#include <stdio.h>      /* printf */
#include <math.h>       /* log */
#include <iostream>
#include <fstream>
#include "mcmc.hpp"




int main () {
    Data data;
    data.calc_mean_predictors();
    Param parameters;
    MCMC mcmc;
    mcmc.run_mcmc(parameters, data);
    return 0;
}
