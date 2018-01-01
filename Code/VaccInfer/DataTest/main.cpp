//
//  ParamTest.cpp
//  VaccInfer
//
//  Created by Lucy Li on 12/21/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#include "../VaccInfer/data.hpp"


int main () {
    int num_ind = 500;
    int num_types = 31;
    int num_times = 4;
    double times_arr[4] = {182, 213, 365, 395};
    std::vector <double> times;
    for (int i=0; i<num_times; i++) times.push_back(times_arr[i]);
    int num_predictors = 13;
    Data data(num_ind, num_types, num_times, times, num_predictors);
    std::cout << "Empty instance of Data generated." << std::endl;
    std::cout << "Number of individuals: " << data.n_ind << std::endl;
    std::cout << "Number of serotypes: " << data.n_tot << std::endl;
    std::cout << "Number of time points: " << data.n_time << std::endl;
    if (data.n_time != data.times.size()) {
        std::cout << "ERROR: n_time does not equal the vector of times" << std::endl;
        std::cout << "n_time = " << data.n_time << std::endl;
        std::cout << "data.times.size() = " << data.times.size() << std::endl;
        return 1;
    }
    if ((data.n_ind * data.n_time)!=data.carriage.size()) {
        std::cout << "ERROR: carriage data is the wrong size. Should be n_ind x n_time" << std::endl;
        return 1;
    }
    std::cout << "Carriage data length: n_ind x n_time = " << data.carriage.size() << std::endl;
    if ((data.n_ind * data.n_predictors)!=data.predictors.size()) {
        std::cout << "ERROR: predictors matrix is the wrong size. Should be n_ind x n_predictors" << std::endl;
        return 1;
    }
    std::cout << "Predictors length: num_ind x num_predictors: " << data.predictors.size() << std::endl;
    if ((data.n_ind * data.n_predictors)!=data.predictor_map.size()) {
        std::cout << "ERROR: predictors map is the wrong size. Should be n_ind x n_predictors" << std::endl;
        return 1;
    }
    std::cout << "Predictors map length: num_ind x num_predictors: " << data.predictor_map.size() << std::endl;
    std::cout << "DATA SUCCESSFULLY GENERATED." << std::endl;
    return 0;
}


