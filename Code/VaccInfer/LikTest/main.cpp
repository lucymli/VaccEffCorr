//
//  main.cpp
//  LikTest
//
//  Created by Lucy Li on 12/27/17.
//  Copyright © 2017 Lucy M Li, CCDD, HPSH. All rights reserved.
//

#include <assert.h>
#include "../VaccInfer/likelihood.hpp"

void test_prediction_func (bool verbose=false) {
    std::vector <double> vals = {-2.5};
    while (vals.back() < 2.5) vals.push_back(vals.back() + 0.5);
    std::vector <double> a = {0};
    while (a.back() < 2) a.push_back(a.back() + 0.2);
    std::vector <double> b = {-2};
    while (b.back() < 0) b.push_back(b.back() + 0.2);
    std::cout << "Testing prediction_func" << std::endl;
    std::vector <bool> matched;
    std::vector <double> computed, expected;
    for (int i=0; i<vals.size(); i++) {
        for (int j=0; j<a.size(); j++) {
            for (int k=0; k<b.size(); k++) {
                computed.push_back(prediction_func(vals[i], a[j], b[k]));
                expected.push_back(std::min(1.0, a[j]*exp(vals[i]*b[k])));
                if (verbose) {
                    std::cout << "val = " << vals[i] << "; a = " << a[j] << "; b = " << b[k] << std::endl;
                    std::cout << "Computed: " << computed.back();
                    std::cout << "\tExpected: " << expected.back() << std::endl;
                } else {
                    matched.push_back(std::abs(computed.back()-expected.back())<1e-8);
                }
            }
        }
    }
    if (std::all_of(matched.begin(), matched.end(), [](bool v) { return v; })) {
        std::cout << "TEST SUCCESS. All predicted values are correct." << std::endl;
    } else {
        std::cout << "TEST FAILURE. The following predictions were incorrect:" << std::endl;
        for (int i=0; i<matched.size(); i++) {
            if (!matched[i]) {
                std::cout << "val = " << vals[i/(a.size()*b.size())] << "; a = " << a[(i%(a.size()*b.size()))/b.size()] << "; b = " << b[i%b.size()] << std::endl;
                std::cout << "Computed: " << computed[i];
                std::cout << "\tExpected: " << expected[i] << std::endl;
            }
        }
    }
    std::cout << "====================================" << std::endl;
}

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
    std::cout << "\tExpected: " << std::endl;
    std::cout << "-0.4, -2.5, -0.4" << std::endl;
    printf("====================================\n");
    //
    // Test the prediction_func
    //
    test_prediction_func ();
    arma::mat base (parameters.n_tot+1, parameters.n_tot+1);
    return 0;
}
