//
//  mcmc.cpp
//  VaccInfer
//
//  Created by Lucy Li on 5/16/17.
//  Copyright © 2017 Lucy M Li, CCDD, HPSH. All rights reserved.
//

#include "mcmc.hpp"


MCMC::MCMC() {
    iter = 0;
    niter = 100000;
    sample_every = 100;
    adapt_optimal=0.23;
    adapt_every = 10;
    adapt_until = 100;
    tot_blocks = 5;
}
