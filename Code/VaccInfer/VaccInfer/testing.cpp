//
//  testing.cpp
//  VaccInfer
//
//  Created by Lucy Li on 12/30/17.
//  Copyright © 2017 Lucy M Li, CCDD, HPSH. All rights reserved.
//

#include "testing.hpp"

void horizontal_rule () {
    std::cout << "====================================" << std::endl;
}

void print_matrix (arma::mat mat_to_output, std::string filename, int mat_dim) {
    std::ofstream outputfile;
    outputfile.open(filename);
    for (int row_i=0; row_i<=mat_dim; row_i++) {
        outputfile << mat_to_output(row_i, 0);
        for (int col_i=1; col_i<=mat_dim; col_i++) {
            outputfile << "\t" << mat_to_output(row_i, col_i);
        }
        outputfile << std::endl;
    }
    outputfile.close();
}

