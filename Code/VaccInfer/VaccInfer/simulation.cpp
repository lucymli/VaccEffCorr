//
//  simulation.cpp
//  VaccInfer
//
//  Created by Lucy Li on 12/21/17.
//  Copyright © 2017 Lucy M Li, CCDD, HPSH. All rights reserved.
//

#include "simulation.hpp"

void Simulation::Simulation(Param parameters, Data sim_data) {
    int nthread = 4;//OMP_NUM_THREADS;
    arma::mat base (parameters.n_tot+1, parameters.n_tot+1);
    fill_rates(parameters, base);
    double time_diff;
    int previous, now;
    #pragma omp parallel for schedule(static, 1)
    for (int tn=0; tn<nthread; tn++) {
        arma::mat exp_mat = arma::expmat(base*sim_data.times[0]);
        std::vector <double> cumulative (1, exp_mat(0, 0));
        for (int i=1; i<=n_tot; i++) cumulative.push_back(cumulative[i-1]+exp_mat(0, i));
        for (int ind_i=tn; ind_i<sim_data.n_ind; ind_i+=nthread) {
            double rnum = runif();
            double x = cumulative[0];
            int event = 0;
            while (rnum < x) {
                event++;
                x = cumulative[event];
                if (event >= (parameters.n_tot + 1)) break;
            }
            sim_data.set_carriage(ind_i, 0, (double) event-1);
        }
    }
    for (int time_t=1; time_t<data.n_time; time_t++) {
#pragma omp parallel for schedule(static, 1)
        for (int tn=0; tn<nthread; tn++) {
            double time_interval = sim_data.times[time_t] - sim_data.times[time_t-1];
            arma::mat exp_mat = arma::expmat(base*time_interval);
            for (int ind_i=tn; ind_i<sim_data.n_ind; ind_i+=nthread) {
                int prev = sim_data.get_carriage(ind_i, time_t-1);
                std::vector <double> cumulative (1, exp_mat(prev, 0));
                for (int i=1; i<=n_tot; i++) cumulative.push_back(cumulative[i-1]+exp_mat(prev, i));
                double rnum = runif();
                double x = cumulative[0];
                int event = 0;
                while (rnum < x) {
                    event++;
                    x = cumulative[event];
                    if (event >= (parameters.n_tot + 1)) break;
                }
                sim_data.set_carriage(ind_i, 0, (double) event-1);
            }
        }
    }
}
