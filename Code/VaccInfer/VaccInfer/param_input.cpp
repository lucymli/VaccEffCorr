//
//  param.cpp
//  VaccInfer
//
//  Created by Lucy Li on 2/13/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//


#include "param.hpp"


Param::Param () {
    SMALLEST_NUMBER = -std::numeric_limits<float>::max()+100000.0;
    output_file_name = "output.txt";
    n_vtypes = 13;
    n_tot = 33;
    n_params = 99;
    n_blocks = 6;
    llik = 0.0;
    new_llik = 0.0;
    lprior = 0.0;
    new_lprior = 0.0;
    block_ptr = 0;
    n_ind = 200;
    inferred_risk_file = "inferred_risk.txt";
    inferred_risk.resize(n_ind*n_vtypes); // n_ind * n_vtypes
    inferred_risk_temp.resize(n_ind*n_vtypes);
    int num_per_block_array[6] = {35, 35, 2, 13, 13, 1};
    std::copy(num_per_block_array, num_per_block_array+n_blocks, std::back_inserter(num_per_block));
    block_starts.push_back(0);
    for (int i=1; i<num_per_block.size(); i++) block_starts.push_back(num_per_block[i-1]+block_starts[i-1]);
    double params_array[99] =
    {
        0.001315481, 0.005084575, // mean and standard deviation of lambda
        0.000535396, 0.00201451, 0.000249904, 0.00153199, 0.000490381, 0.00115203, 0.00124178, 0.000977398, 0.000123424, 0.000314247, 0.000082313, 0.000177742, 0.000130397,
        0.0000137666, 0.00000189964, 0.000147169, 0.000834131, 0.0010559, 0.000693417, 0.000428734, 0.0000244788, 0.000189081, 0.000508208, 0.00000119155, 0.0000903176, 0.000118967, 0.0000688422, 0.0000459625, 0.0000292863, 0.000625802, 0.00000502775, 0.00000806392, 0.0294991,
        0.02381501, 0.01531154, // mean and standard deviation of mu
        0.0429771, 0.000417408, 0.00170975, 0.0374325, 0.0108858, 0.000289974, 0.0137252, 0.0392182, 0.000231978, 0.0362689, 0.00627818, 0.0348528, 0.0494624,
        0.0176524, 0.042872, 0.0458283, 0.01483, 0.0137507, 0.0340181, 0.0386283, 0.025324, 0.0167337, 0.0148374, 0.0378644, 0.0321509, 0.0436316, 0.0206733, 0.00467575, 0.0353994, 0.0207082, 0.0192338, 0.0273504, 0.00598254,
        0.429731, 0.379037, // competition parameters
        0.00493, 0.04028, 0.01644, 0.0537, 0.01863, 0.08408, 0.04125, 0.0228, 0.01228, 0.02089, 0.01839, 0.00417, 0.00848, // first parameter of antibody model
        -0.00648, -0.13289, 3.98598, -0.02095, -0.06409, 0.29593, 0.16986, -0.05213, -0.10403, 0.60803, 1.12495, -0.00274, -0.03166, // second parameter of antibody model
        1 // surrogate indicator
    };
    for (int i=0; i<n_params; i++) {
        params.push_back(params_array[i]);
        tempparam.push_back(params_array[i]);
    }
    double params_sd_array[6] = {0.00005, 0.0001, 0.05, 0.005, 0.05, 1};
    for (int i=0; i<n_blocks; i++) params_sd.push_back(params_sd_array[i]);
    accepted.resize(n_blocks, 0);
    rejected.resize(n_blocks, 0);
}

double Param::calc_lprior() {
    double dens = 0.0;
    // lambda
//    dens += dnorm(params[0], 0.001315481, 0.0013);
//    dens += dnorm(params[1], 0.005084575, 0.005084575);
    for (int i=block_starts[0]; i<block_starts[0]+num_per_block[0]; i++) dens += dunif(params[i], 0.000001, 0.1);
    // mu
//    dens += dnorm(params[0], 0.02381501, 0.02381501);
//    dens += dnorm(params[1], 0.01531154, 0.01531154);
    for (int i=block_starts[1]; i<block_starts[1]+num_per_block[1]; i++) dens += dunif(params[i], 0.000001, 0.1);
    // competition
    for (int i=block_starts[2]; i<block_starts[2]+num_per_block[2]; i++) dens += dunif(params[i], 0.0, 1.0);
    // antibody model parameter 1
    for (int i=block_starts[3]; i<block_starts[3]+num_per_block[3]; i++) dens += dunif(params[i], 0.0, 1.0);
    // antibody model parameter 2
    for (int i=block_starts[4]; i<block_starts[4]+num_per_block[4]; i++) dens += dunif(params[i], -5.0, 5.0);
    // surrogate indicator for 6B
//    if (params[block_starts[5]] < 0 | params[block_starts[5]] >= n_vtypes) dens += SMALLEST_NUMBER;
    if (dens < SMALLEST_NUMBER) dens = SMALLEST_NUMBER;
    return (dens);
}

void Param::propose() {
    int ntries = 50;
    tempparam.swap(params);
    if (block_ptr < 2) {
//        params[block_starts[0]] = rnorm(tempparam[block_starts[0]], params_sd[block_ptr], 0.0, 1.0, ntries);
//        params[block_starts[0]+1] = rnorm(tempparam[block_starts[1]], params_sd[block_ptr], 0.0, 1.0, ntries);
        for (int i=block_starts[block_ptr]; i<block_starts[block_ptr]+num_per_block[block_ptr]; i++) {
//            params[i] = rnorm(params[block_starts[block_ptr]], params[block_starts[block_ptr]+1], 0.0, 1.0, ntries);
              params[i] = rnorm(tempparam[i], params_sd[block_ptr], 0.0, 1.0, ntries);

        }
    }
    if (block_ptr==2) {
        for (int i=block_starts[block_ptr]; i<block_starts[block_ptr]+num_per_block[block_ptr]; i++) {
            params[i] = rnorm(tempparam[i], params_sd[block_ptr], 0.0, 1.0, ntries);
        }
    }
    if (block_ptr==3) {
        for (int i=block_starts[block_ptr]; i<block_starts[block_ptr]+num_per_block[block_ptr]; i++) {
            params[i] = rnorm(tempparam[i], params_sd[block_ptr], 0.0, 1.0, ntries);
        }
    }
    if (block_ptr==4) {
        for (int i=block_starts[block_ptr]; i<block_starts[block_ptr]+num_per_block[block_ptr]; i++) {
            params[i] = rnorm(tempparam[i], params_sd[block_ptr], -5.0, 5.0, ntries);
        }
    }
    if (block_ptr==5) {
        params[block_starts[block_ptr]] = runif(0.0, (double) n_vtypes);
    }
    inferred_risk.swap(inferred_risk_temp);
}

void Param::initialize_file () {
    std::ofstream output_file;
    output_file.open(output_file_name);
    output_file << "state\tposterior\tlikelihood\tprior";
    for (int i=0; i<n_tot; i++) output_file << "\tlambda" << i;
    for (int i=0; i<n_tot; i++) output_file << "\tmu" << i;
    for (int i=0; i<2; i++) output_file << "\tcompetition" << i;
    for (int i=0; i<n_vtypes; i++) output_file << "\tab1type" << i;
    for (int i=0; i<n_vtypes; i++) output_file << "\tab2type" << i;
//    output_file << "\tsurrogate";
    output_file << std::endl;
    output_file.close();
}

void Param::write_inferred_risk (int iter) const {
    std::ofstream outputfile;
    if (iter==0) {
        outputfile.open(inferred_risk_file);
        outputfile << "state\t";
        for (int ind_i=0; ind_i<n_ind; ind_i++) {
            for (int type_i=0; type_i<n_vtypes; type_i++) {
                outputfile << "ind" << ind_i << "type" << type_i << "\t";
            }
        }
        outputfile << std::endl;
    }
    else {
        outputfile.open(inferred_risk_file, std::ios::app);
    }
    outputfile << iter << "\t";
    for (int ind_i=0; ind_i<n_ind; ind_i++) {
        for (int type_i=0; type_i<n_vtypes; type_i++) {
            outputfile << get_inferred_risk(ind_i, type_i) << "\t";
        }
    }
    outputfile << std::endl;
    outputfile.close();
}
