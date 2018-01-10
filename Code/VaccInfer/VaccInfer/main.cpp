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

//#include <cstdio>      /* printf */
#include <math.h>       /* log */
#include <iostream>
#include <fstream>
#include "../VaccInfer/mcmc.hpp"




int main (int argc, char *argv[]) {
    std::string dataf (argv[2]), paramf(argv[3]), mcmcf(argv[4]);
    Data data(dataf);
//    data.calc_mean_predictors();
    Param parameters(paramf);
    MCMC mcmc(mcmcf);
    mcmc.run_mcmc(parameters, data);
    return 0;
}
