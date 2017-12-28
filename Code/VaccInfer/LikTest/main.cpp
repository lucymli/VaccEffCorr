//
//  main.cpp
//  LikTest
//
//  Created by Lucy Li on 12/27/17.
//  Copyright © 2017 Lucy M Li, CCDD, HPSH. All rights reserved.
//

#include <assert.h>
#include "../VaccInfer/likelihood.hpp"


int main(int argc, const char * argv[]) {
    int num_ind = 500;
    int num_types = 31;
    int num_times = 4;
    double times_arr[4] = {182, 213, 365, 395};
    std::vector <double> times;
    for (int i=0; i<num_times; i++) times.push_back(times_arr[i]);
    int num_predictors = 13;
    Data data(num_ind, num_types, num_times, times, num_predictors);
    Param parameters("../ParamTest/test_params.txt");
    //
    // Test the set_diag_as_negrowsum function
    //
    arma::mat simple_mat = {
        {0.5, 0.3, 0.1},
        {2.5, 0.3, -1.0},
        {0.4, 0.0, 0.9}
    };
    set_diag_as_negrowsum(simple_mat);
    std::cout << "Testing set_diag_as_negrowsum " << std::endl;
    std::cout << "Computed: " << std::endl;
    std::cout << simple_mat(0,0) << ", " << simple_mat(1,1) << ", " << simple_mat(2,2) << std::endl;
    std::cout << "Expected: " << std::endl;
    std::cout << "-0.4, -2.5, -0.4" << std::endl;
    printf("====================================\n");
    //
    // Test the prediction_func
    //
    std::cout << "Testing prediction_func" << std::endl;
    double val = -2.0, a = 2.0, b = -2.0;
    std::cout << "val = " << val << "; a = " << a << "; b = " << b << std::endl;
    std::cout << "Computed: " << prediction_func(val, a, b) << std::endl;
    std::cout << "Expected: " << std::min(1.0, a*exp(val*b)) << std::endl;
    arma::mat base (parameters.n_tot+1, parameters.n_tot+1);
    return 0;
}
