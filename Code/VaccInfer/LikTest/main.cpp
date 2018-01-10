//
//  main.cpp
//  LikTest
//
//  Created by Lucy Li on 12/27/17.
//  Copyright © 2017 Lucy M Li, CCDD, HPSH. All rights reserved.
//

#include "../VaccInfer/simulation.hpp"
#include "../VaccInfer/likelihood.hpp"
#include "../VaccInfer/testing.hpp"

bool check_all_rows_sum_to_zero (arma::mat matrix, bool verbose) {
    set_diag_as_negrowsum(matrix);
    arma::colvec sumcheck = sum(matrix, 1);
    bool notzero = false;
    for (int x=0; x<sumcheck.size(); x++) {
        notzero = sumcheck[x] > 1e-8;
        if (notzero) {
            std::cout << "Row " << x << " does not sum to zero." << std::endl;
            break;
        }
    }
    if (verbose) {
        if (!notzero) std::cout << "TEST SUCCESS. All rows sum to 0." << std::endl;
    }
    return !notzero;
}

void test_set_diag_as_negrowsum () {
    int nrow;
    int num_tests = 100;
    std::vector <bool> test_results;
    std::cout << "Testing set_diag_as_negrowsum" << std::endl;
    for (int i=0; i<num_tests; i++) {
        arma::mat simple_mat;
        if (i<=1) {
            nrow = 3;
            simple_mat = arma::zeros(nrow, nrow);
            simple_mat.fill(i);
        }
        else {
            nrow = rsample(3, 200);
            simple_mat = arma::randu(nrow, nrow);
        }
        test_results.push_back(check_all_rows_sum_to_zero(simple_mat, false));
    }
    horizontal_rule();
}

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
    horizontal_rule();
}

int main(int argc, const char * argv[]) {
    int num_ind = 500;
    int num_types = 31;
    int num_times = 4;
    std::vector <double> times = {182, 213, 365, 395};
    int num_predictors = 13;
    Param parameters("../ParamTest/test_params.txt");
    //
    // Test the set_diag_as_negrowsum function
    //
    test_set_diag_as_negrowsum();
    //
    // Test the prediction_func
    //
    test_prediction_func ();
    //
    // Test fill_rates
    //
    arma::mat base (parameters.n_tot+1, parameters.n_tot+1);
    parameters.fill_rates(base);
    std::cout << "Testing fill_rates:" << std::endl;
    bool check_result = check_all_rows_sum_to_zero (base, true);
    horizontal_rule();
    print_matrix(base, "LikTestFillRatesMatrix.txt", parameters.n_tot);
    //
    // Test calc_llik()
    //
    Data sim_data(num_ind, num_types, num_times, times, num_predictors);
    simulate (parameters, sim_data, true);
    double llik = calc_llik(parameters, sim_data, false);
    std::cout << "Testing calc_llik:" << std::endl;
    std::cout << "Log Likelihood = " << llik << std::endl;
    return 0;
}
