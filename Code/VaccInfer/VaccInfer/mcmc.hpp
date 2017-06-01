//
//  mcmc.hpp
//  VaccInfer
//
//  Created by Lucy Li on 5/16/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#ifndef mcmc_hpp
#define mcmc_hpp

#include <stdio.h>
#include "param.hpp"
#include "likelihood.hpp"

class MCMC {
public:
    MCMC ();
    int iter;
    int niter;
    int sample_every;
    double adapt_optimal;
    int adapt_every;
    int adapt_until;
    int tot_blocks;
    void run_mcmc(Param, Data);
};

#endif /* mcmc_hpp */
