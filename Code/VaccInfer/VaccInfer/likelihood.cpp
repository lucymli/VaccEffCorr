//
//  likelihood.cpp
//  VaccInfer
//
//  Created by Lucy Li on 5/8/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#include "likelihood.hpp"







double calc_llik (Param parameters_input, Data data_input, bool use_mean_ab) {
    int nthread = 4;//OMP_NUM_THREADS;
    //#pragma omp parallel for schedule(static, 1)
    std::vector <Param> param_vec;
    std::vector <Data> data_vec;
    std::vector <double> llik_vec;
    for (int i=0; i<nthread; i++) {
        param_vec.push_back(parameters_input);
        data_vec.push_back(data_input);
        llik_vec.push_back(0.0);
    }
    for (int tn=0; tn<nthread; tn++) {
        double lik;
        Param parameters = param_vec[tn];
        Data data = data_vec[tn];
        arma::mat base (parameters.n_tot+1, parameters.n_tot+1);
        parameters.fill_rates(base);
        double time_diff;
        int previous, now;
        for (int ind_i=tn; ind_i<data.n_ind; ind_i+=nthread) {
            arma::mat ind_mat (base);
            for (int time_t=0; time_t<data.n_time; time_t++) {
                if (time_t > 0) {
                    parameters.predict_lambda(ind_mat, data, ind_i, use_mean_ab);
                    time_diff = data.times[time_t]-data.times[time_t-1];
                    previous = data.get_carriage(ind_i, time_t-1);
                }
                else {
                    time_diff = data.times[time_t];
                    previous = 0;
                }
                now = data.get_carriage(ind_i, time_t);
                arma::mat exp_mat = arma::expmat(ind_mat*time_diff);
                lik = exp_mat(previous, now);
                llik_vec[tn] += std::log10(lik);
            }
        }
    }
    double llik = std::accumulate(llik_vec.begin(), llik_vec.end(), 0.0);
    if ((llik < parameters_input.SMALLEST_NUMBER) | (!std::isfinite(llik))) llik = parameters_input.SMALLEST_NUMBER;
    return (llik);
}

