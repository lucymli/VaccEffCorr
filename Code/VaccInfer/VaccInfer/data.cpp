//
//  data.cpp
//  VaccInfer
//
//  Created by Lucy Li on 2/13/17.
//  Copyright © 2017 Lucy M Li, CCDD, HPSH. All rights reserved.
//

#include "data.hpp"


double Data::get_carriage(int ind_i, int time_i) const {
    return (carriage[time_i * n_ind + ind_i]);
}

double Data::get_metadata(int ind_i, int i) const {
    return (metadata[i*n_ind + ind_i]);
}

double Data::get_predictor(int ind_i, int predict_i) const {
    return (predictors[predict_i * n_ind + ind_i]);
}

double Data::get_inferred_risk(int ind_i, int i) const {
    return (inferred_risk[i * n_ind + ind_i]);
}

double Data::sum_risk_by_ind(int ind_i) const {
    double risk=0.0;
    for (int i=0; i<n_vacc_types; i++) {
        risk += inferred_risk[i*n_ind + ind_i];
    }
    return (risk);
}

double Data::sum_risk_by_type(int type_i) const {
    double risk=0.0;
    for (int i=0; i<n_ind; i++) {
        risk += inferred_risk[type_i*n_ind+i];
    }
    return (risk);
}

double Data::mean_risk_by_ind(int ind_i) const {
    double mean_risk = this->sum_risk_by_ind(ind_i);
    mean_risk /= (double) n_vacc_types;
    return (mean_risk);
}

double Data::mean_risk_by_type(int type_i) const {
    double mean_risk = this->sum_risk_by_type(type_i);
    mean_risk /= (double) n_ind;
    return (mean_risk);
}

