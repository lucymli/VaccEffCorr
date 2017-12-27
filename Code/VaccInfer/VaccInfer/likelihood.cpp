//
//  likelihood.cpp
//  VaccInfer
//
//  Created by Lucy Li on 5/8/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#include "likelihood.hpp"

void set_diag_as_negrowsum (arma::mat & rates) {
    for (int i=0; i<rates.n_rows; i++) {
        rates(i, i) =  - (accu(rates.row(i)) - rates(i, i));
    }
}

double prediction_func (double val, double a, double b) {
    double result = a * exp(val * b);
    if (result < 0) result = 0.0;
    if (result > 1) result = 1.0;
    return result;
}

void predict_lambda (arma::mat & rates, Param &parameters, Data &data, int ind_i, bool use_mean_ab) {
    double a, b;
    double predictor;
    int pos;
    double multiplier;
    for (int type_i=0; type_i<parameters.n_vtypes; type_i++) {
        a = parameters.params[parameters.n_tot*2+type_i];
        b = parameters.params[parameters.n_tot*2+parameters.n_vtypes+type_i];
        pos = data.get_predictor_index(ind_i, type_i);
        if (pos < 0) predictor = 1.0;
        else {
            if (use_mean_ab) predictor = data.mean_predictors[ind_i];
            else {
                predictor = data.get_predictor(ind_i, pos);
            }
            multiplier = prediction_func(predictor, a, b);
            if (multiplier < 1.0) {
                rates.col(type_i+1) *= multiplier;
            }
            set_diag_as_negrowsum(rates);
        }
    }
    
}

void fill_rates (Param parameters, arma::mat & mat) {
    for (int i=0; i<parameters.n_tot; i++) {
        mat(i+1, 0) = parameters.params[parameters.n_tot+i]; // mu
        mat(0, i+1) = parameters.params[i]; // lambda
    }
    double competition;
    for (int row_i=1; row_i<=parameters.n_tot; row_i++) {
        for (int col_i=1; col_i<=parameters.n_tot; col_i++) {
            if (row_i!=col_i) {
                if (col_i <= parameters.n_vtypes) competition = parameters.params[parameters.n_tot*2+col_i-1];
                else competition = parameters.params[parameters.n_tot*2+parameters.n_vtypes];
                mat(row_i, col_i) = parameters.params[col_i-1]*competition;
            }
        }
    }
    set_diag_as_negrowsum(mat);
}


double calc_llik (Param &parameters, Data &data, bool use_mean_ab) {
    int nthread = 4;//OMP_NUM_THREADS;
    std::vector <double> llik_vec(nthread, 0.0);
    arma::mat base (parameters.n_tot+1, parameters.n_tot+1);
    fill_rates(parameters, base);
    double time_diff;
    int previous, now;
    //#pragma omp parallel for schedule(static, 1)
    for (int tn=0; tn<nthread; tn++) {
        for (int ind_i=tn; ind_i<data.n_ind; ind_i+=nthread) {
            predict_lambda(parameters, data, use_mean_ab);
            for (int time_t=1; time_t<data.n_time; time_t++) {
                time_diff = data.times[time_t]-data.times[time_t-1];
                arma::mat ind_mat = arma::expmat(fill_rates(base, data, parameters, ind_i)*time_diff);
                previous = data.get_carriage(ind_i, time_t-1);
                now = data.get_carriage(ind_i, time_t);
                llik_vec[tn] += std::log10(ind_mat[previous, now]);
            }
        }
    }
    double llik = std::accumulate(llik_vec.begin(), llik_vec.end(), 0.0);
    if ((llik < parameters.SMALLEST_NUMBER) | (!std::isfinite(llik))) llik = parameters.SMALLEST_NUMBER;
    return (llik);
}
