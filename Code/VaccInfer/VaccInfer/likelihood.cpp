//
//  likelihood.cpp
//  VaccInfer
//
//  Created by Lucy Li on 5/8/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#include "likelihood.hpp"




double predict_lambda (Param &parameters, Data data, bool use_mean_ab) {
    double a, b;
    double predictor;
    for (int type_i=0; type_i<parameters.n_vtypes; type_i++) {
        a = parameters.params[parameters.block_starts[3]+type_i];
        b = parameters.params[parameters.block_starts[4]+type_i];
        for (int ind_i=0; ind_i<data.n_ind; ind_i++) {
            if (use_mean_ab) predictor = data.mean_predictors[ind_i];
            else predictor = data.get_predictor(ind_i, type_i);
            parameters.set_inferred_risk(a*exp(b*predictor), ind_i, type_i);//Exp
//            parameters.set_inferred_risk(a+b*predictor, ind_i, type_i);//Linear
        }
    }
}

void fill_rates (Param parameters, arma::mat & mat) {
    double tot_rate = 0.0;
    for (int i=0; i<parameters.n_tot; i++) {
        mat(i+1, 0) = parameters.get_param(1, i); // mu
        mat(0, i+1) = parameters.get_param(0, i); // lambda
        tot_rate += mat(0, i+1);
    }
    double competition;
    for (int row_i=1; row_i<=parameters.n_tot; row_i++) {
        tot_rate = mat(row_i, 0);
        for (int col_i=1; col_i<=parameters.n_tot; col_i++) {
            if (row_i!=col_i) {
                if (col_i <= parameters.n_vtypes) competition = parameters.get_param(2, col_i-1);
                else competition = parameters.get_param(2, parameters.n_vtypes);
                mat(row_i, col_i) = parameters.get_param(0, col_i-1)*competition;
                tot_rate += mat(row_i, col_i);
            }
        }
        mat(row_i, row_i) = -tot_rate;
    }
}

arma::mat fill_rates (arma::mat base_mat, Data data, Param parameters, int ind_i) {
    arma::mat new_mat(base_mat);
    for (int row_i=0; row_i<=parameters.n_vtypes; row_i++) {
        for (int col_i=1; col_i<=parameters.n_vtypes; col_i++) {
            if (row_i!=col_i) {
//                if (col_i==2) {
//                    int choice_of_correlate = (int)parameters.get_param(5, 0);
//                    new_mat(row_i, col_i) *= parameters.get_inferred_risk(ind_i, choice_of_correlate);
//                }
//                else new_mat(row_i, col_i) *= parameters.get_inferred_risk(ind_i, col_i-1);
                if (ind_i < 642) {
                    if (col_i < 7) new_mat(row_i, col_i) *= parameters.get_inferred_risk(ind_i, col_i-1);
                    if (col_i == 11) new_mat(row_i, col_i) *= parameters.get_inferred_risk(ind_i, 2);
                }
                else {
                    new_mat(row_i, col_i) *= parameters.get_inferred_risk(ind_i, col_i-1);
                }
            }
        }
    }
    double tot_rate;
    for (int row_i=0; row_i<=parameters.n_tot; row_i++) {
        tot_rate = 0.0;
        for (int col_i=0; col_i<parameters.n_tot; col_i++) {
            if (row_i!=col_i) tot_rate += new_mat(row_i, col_i);
        }
        new_mat(row_i, row_i) = -tot_rate;
    }
    return (new_mat);
}


double calc_llik (Param &parameters, Data &data, bool use_mean_ab) {
    int nthread = 4;//OMP_NUM_THREADS;
    predict_lambda(parameters, data, use_mean_ab);
    std::vector <double> llik_vec(nthread, 0.0);
    arma::mat base (parameters.n_tot+1, parameters.n_tot+1);
    fill_rates(parameters, base);
    double time_diff;
    int previous, now;
    for (int time_t=1; time_t<data.n_time; time_t++) {
        #pragma omp parallel for schedule(static, 1)
        for (int tn=0; tn<nthread; tn++) {
            for (int ind_i=tn; ind_i<data.n_ind; ind_i+=nthread) {
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
