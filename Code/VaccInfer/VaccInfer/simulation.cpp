//
//  simulation.cpp
//  VaccInfer
//
//  Created by Lucy Li on 12/21/17.
//  Copyright © 2017 Lucy M Li, CCDD, HPSH. All rights reserved.
//

#include "simulation.hpp"

void simulate_predictor (Param & parameters, Data &sim_data, std::vector <double> predictor_means, std::vector <double> predictor_sd) {
    double predictor_value;
    for (int ind_i=0; ind_i<sim_data.n_ind; ind_i++) {
        for (int predictor_i=0; predictor_i<sim_data.n_predictors; predictor_i++) {
            predictor_value = rnorm(predictor_means[predictor_i], predictor_sd[predictor_i]);
            sim_data.set_predictor(ind_i, predictor_i, predictor_value);
            sim_data.set_predictor_index(ind_i, predictor_i, predictor_i);
        }
    }
}

void simulate (Param &parameters_input, Data &sim_data_input, bool use_predictors) {
    int nthread = 4;//omp_get_max_threads() ;
    arma::mat base (parameters_input.n_tot+1, parameters_input.n_tot+1);
    parameters_input.fill_rates(base);
    std::vector <Param> parameters_vec;
    std::vector <Data> data_vec;
    for (int i=0; i<nthread; i++) {
        parameters_vec.push_back(parameters_input);
        data_vec.push_back(sim_data_input);
    }
    //#pragma omp parallel for schedule(static, 1)
    for (int tn=0; tn<nthread; tn++) {
        Param parameters = parameters_vec[tn];
        Data * sim_data = &data_vec[tn];
        double time_interval;
        int prev;
        for (int ind_i=tn; ind_i<sim_data->n_ind; ind_i+=nthread) {
            arma::mat ind_mat (base);
            if (use_predictors) parameters.predict_lambda(ind_mat, *sim_data, ind_i, false);
            for (int time_t=0; time_t<sim_data->n_time; time_t++) {
                if (time_t>0) {
                    time_interval = sim_data->times[time_t] - sim_data->times[time_t-1];
                    prev = sim_data->get_carriage(ind_i, time_t-1);
                }
                else {
                    time_interval = sim_data->times[time_t];
                    prev = 0;
                }
                arma::mat exp_mat = arma::expmat(base*time_interval);
                std::vector <double> cumulative (1, exp_mat(prev, 0));
                for (int i=1; i<=sim_data->n_tot; i++) cumulative.push_back(cumulative[i-1]+exp_mat(prev, i));
                double rnum = runif();
                double x = cumulative[0];
                int event = 0;
                while (x < rnum) {
                    event++;
                    x = cumulative[event];
                    if (event >= parameters.n_tot) break;
                }
                sim_data->set_carriage(ind_i, time_t, (double) event);
            }
        }
    }
    for (int tn=0; tn<nthread; tn++) {
        Data * sim_data = &data_vec[tn];
        for (int ind_i=tn; ind_i<sim_data->n_ind; ind_i+=nthread) {
            for (int time_t=0; time_t<sim_data->n_time; time_t++) {
                sim_data_input.set_carriage(ind_i, time_t, sim_data->get_carriage(ind_i, time_t));
            }
        }
    }
}
