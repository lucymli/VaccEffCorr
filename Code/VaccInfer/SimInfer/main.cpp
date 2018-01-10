//
//  main.cpp
//  SimInfer
//
//  Created by Lucy Li on 1/10/18.
//  Copyright © 2018 Lucy M Li, CCDD, HPSH. All rights reserved.
//

#include <iostream>
#include "../VaccInfer/simulation.hpp"


int main(int argc, const char * argv[]) {
    // args: 1) -I.. 2) param.txt 3) data.txt
    std::string input_param_file (argv[2]);
    Param parameters(input_param_file);
    std::vector <double> times = {182, 213, 365, 395, 548, 730};
    Data sim_data(parameters.n_ind, parameters.n_tot, 6, times, parameters.n_vtypes);
    simulate(parameters, sim_data, false);
    std::vector <double> predictor_means = {0.337897650192541, 0.358080135848048, 0.147576770946833, 0.765048005860137, 0.178581563922821, 0.463256490662169, 0.0593742682792419, 0.320398220824545, -0.0123828538755206, 0.145830330191984, 0.406638314189488, 0.524244394092351, 0.263846099885817};
    std::vector <double> predictor_sd = {0.333371215952248, 0.509179146638132, 0.321682578656864, 0.505434602688631, 0.3685758533035, 0.365853887982019, 0.467887033893633, 0.351904953410131, 0.36599584083075, 0.351151567991839, 0.429909039874033, 0.309507219772883, 0.347467866724146};
    simulate_predictor(parameters, sim_data, predictor_means, predictor_sd);
    simulate(parameters, sim_data, true);
    std::string data_file_name (argv[3]);
    sim_data.print_data_to_file(data_file_name);
    return 0;
}
