//
//  data.cpp
//  VaccInfer
//
//  Created by Lucy Li on 2/13/17.
//  Copyright © 2017 Lucy M Li, CCDD, HPSH. All rights reserved.
//

#include "data.hpp"



Data::Data (std::vector<double> vdata, std::vector<double>nvdata,
            std::vector<double> abdata, std::vector <double>timings,
            int vtypes, int nvtypes, int nvacc, int nnvacc) {
    for (auto i=vdata.begin(); i<vdata.end(); i++) swab_data_v.push_back(*i);
    for (auto i=nvdata.begin(); i<nvdata.end(); i++) swab_data_nv.push_back(*i);
    for (auto i=abdata.begin(); i<abdata.end(); i++) ab_data.push_back(*i);
    swab_times.push_back(timings[0]);
    for (int i=1; i<timings.size(); i++) swab_times.push_back(timings[i]-timings[i-1]);
    n_swabs = timings.size();
    n_vacc = nvacc;
    n_nvacc = nnvacc;
    n_vtypes = vtypes;
    n_nvtypes = nvtypes;
    n_tot = n_vtypes + n_nvtypes;
    for (int i=0; i<n_vtypes; i++) {
        double max_val = abdata[i*n_vacc];
        double min_val = abdata[i*n_vacc];
        double mean_val = max_val;
        for (int ind_i=1; ind_i<n_nvacc; ind_i++) {
            double curr_val = abdata[i*n_vacc+ind_i];
            mean_val += curr_val;
            if (curr_val > max_val) max_val = curr_val;
            if (curr_val < min_val) min_val = curr_val;
        }
        max_ab.push_back(max_val);
        min_ab.push_back(min_val);
        mean_ab.push_back(mean_val/(double)n_vacc);
    }
}

int Data::operator[] (std::string name) const {
    if (name=="n_vacc") return (n_vacc);
    else if (name=="n_nvacc") return (n_nvacc);
    else if (name=="n_swabs") return (n_swabs);
    else if (name=="n_vtypes") return (n_vtypes);
    else if (name=="n_nvtypes") return (n_nvtypes);
    else return(n_tot);
}

double Data::get_swab_time (int time) {
    return (swab_times[time]);
}

double Data::get_swab_v (int ind_i, int swab_i) {
    return (swab_data_v[swab_i*n_vacc+ind_i]);
}
double Data::get_swab_nv (int ind_i, int swab_i) {
    return (swab_data_nv[swab_i*n_nvacc+ind_i]);
}

double Data::get_ab (int ind_i, int type_i) {
    return(ab_data[type_i*n_vacc+ind_i]);
}

double Data::get_mean_ab (int type_i) {
    return(mean_ab[type_i]);
}

double Data::get_max_ab (int type_i) {
    return(max_ab[type_i]);
}

double Data::get_min_ab (int type_i) {
    return(min_ab[type_i]);
}

