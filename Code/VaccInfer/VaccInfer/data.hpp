//
//  data.hpp
//  VaccInfer
//
//  Created by Lucy Li on 2/13/17.
//  Copyright © 2017 Lucy M Li, CCDD, HPSH. All rights reserved.
//

#ifndef data_hpp
#define data_hpp

#include <stdio.h>
#include <vector>
#include <string>

class Data {
    std::vector <double> swab_data_v; // nrow=total number of vaccinated, ncol=swabs
    std::vector <double> swab_data_nv; // nrow=total number of vaccinated, ncol=swabs
    std::vector <double> swab_times; // timing of swabs
    std::vector <double> ab_data; // nrow=total number of vaccinated, ncol=serotypes in a vaccine
    std::vector <double> lambda_reduction; // nrow=total number of vaccinated, ncol=serotypes in a vaccine
    int n_vacc, n_nvacc, n_swabs, n_vtypes, n_nvtypes,n_tot;
    std::vector<double> mean_ab, max_ab, min_ab;
public:
    Data (std::vector<double>, std::vector<double>, std::vector<double>,
          std::vector <double>, int, int, int, int);
    int operator[] (std::string) const;
    double get_swab_time (int);
    double get_swab_v (int, int);
    double get_swab_nv (int, int);
    double get_ab (int, int);
    double get_mean_ab (int);
    double get_max_ab (int);
    double get_min_ab (int);
};

#endif /* data_hpp */
