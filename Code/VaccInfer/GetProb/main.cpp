//
//  main.cpp
//  GetProb
//
//  Created by Lucy Li on 2/1/18.
//  Copyright © 2018 Lucy M Li, CCDD, HPSH. All rights reserved.
//

#include <iostream>
#include "../VaccInfer/testing.hpp"
#include "../VaccInfer/likelihood.hpp"

int main(int argc, const char * argv[]) {
    // args: 1) -I.. 2) param.txt 3) data.txt 4) time_interval 5) bool: assume all susceptible at start 6) output file name
    std::string input_param_file (argv[2]);
    Param parameters(input_param_file);
    std::string input_data_file (argv[3]);
    Data data(input_data_file);
    double time_interval = std::stod(argv[4]);
    bool all_susc = std::stoi(argv[5]);
    arma::mat output (data.n_ind, data.n_tot+1);
    get_multinom_prob(parameters, data, time_interval, output, all_susc);
    std::string output_file (argv[6]);
    print_matrix (output, output_file, data.n_ind, data.n_tot+1);
    return 0;
}
