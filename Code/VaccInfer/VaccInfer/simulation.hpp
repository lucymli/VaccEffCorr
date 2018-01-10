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
#include "param.hpp"

void simulate_predictor (Param &, Data &, std::vector <double>, std::vector <double>);

void simulate (Param &, Data &, bool);


#endif /* simulation_hpp */

