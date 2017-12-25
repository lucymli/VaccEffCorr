//
//  main.cpp
//  SimTest
//
//  Created by Lucy Li on 12/21/17.
//  Copyright © 2017 Lucy M Li, CCDD, HPSH. All rights reserved.
//

#include <iostream>
#include "../VaccInfer/simulation.hpp"

int main(int argc, const char * argv[]) {
    int num_ind = 500;
    int num_types = 31;
    int num_times = 4;
    double times_arr[4] = {182, 213, 365, 395};
    std::vector <double> times;
    for (int i=0; i<num_times; i++) times.push_back(times_arr[i]);
    int num_predictors = 13;
    Data data(num_ind, num_types, num_times, times, num_predictors);
    Param parameters("test_params.txt");
    Simulation sim (Param parameters, Data sim_data);
    return 0;
}
