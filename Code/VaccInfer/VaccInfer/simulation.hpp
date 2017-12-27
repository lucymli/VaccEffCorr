//
//  simulation.hpp
//  VaccInfer
//
//  Created by Lucy Li on 12/21/17.
//  Copyright © 2017 Lucy M Li, CCDD, HPSH. All rights reserved.
//

#ifndef simulation_hpp
#define simulation_hpp

#include <stdio.h>
#include <armadillo>
#include "likelihood.hpp"

class Simulation {
public:
    Simulation ();
    Simulation (Param &, Data &);
};


#endif /* simulation_hpp */

