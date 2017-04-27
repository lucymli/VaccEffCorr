//
//  param.cpp
//  VaccInfer
//
//  Created by Lucy Li on 2/13/17.
//  Copyright © 2017 Lucy M Li, CCDD, HPSH. All rights reserved.
//


#include "param.hpp"

double SMALLEST_NUMBER = -std::numeric_limits<float>::max()+100000.0;
double STATIONARY_TIME = 300.0;
boost::mt19937 rng(NULL);


void print_matrix (arma::mat & mat_to_output, std::string filename, int mat_dim) {
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


Param::Param (int ntypes, int ntot, int nparams, std::vector <double> inputs,
              std::vector <double> input_var, bool mean_ab) {
    n_vtypes = ntypes;
    n_tot = ntot;
    n_nvtypes = n_tot-n_vtypes;
    n_params = nparams;
    tempparam = std::vector <double> (n_tot);
//    frailtySI = inputs[n_tot*2+n_vtypes*3];
//    frailtyIS = inputs[n_tot*2+n_vtypes*3+1];
    for (int i=n_tot*2+n_vtypes*2; i!=nparams; i++) {
        interaction.push_back(inputs[i]);
    }
    for (int i=0; i<n_vtypes; i++) {
        lambda.push_back(inputs[i]);
        mu.push_back(inputs[i+n_tot]);
        thetaSI.push_back(inputs[i+n_tot*2]);
//        thetaIS.push_back(inputs[i+n_tot*2+n_vtypes]);
        p0.push_back(inputs[i+n_tot*2+n_vtypes]);
    }
    for (int i=n_vtypes; i<n_tot; i++) {
        lambda.push_back(inputs[i]);
        mu.push_back(inputs[i+n_tot]);
    }
    /* Check the parameters have been read in correctly */
    // for (auto i=lambda.begin(); i<lambda.end(); i++) Rcout << "lambda " << *i << std::endl;
    // for (auto i=mu.begin(); i<mu.end(); i++) Rcout << "mu " << *i << std::endl;
    // for (auto i=p0.begin(); i<p0.end(); i++) Rcout << "p0 " << *i << std::endl;
    // for (auto i=thetaSI.begin(); i<thetaSI.end(); i++) Rcout << "thetaSI " << *i << std::endl;
    // for (auto i=thetaIS.begin(); i<thetaIS.end(); i++) Rcout << "thetaIS " << *i << std::endl;
    // Rcout << "frailtySI: " << frailtySI << " frailtyIS: " << frailtyIS << " interaction: " << interaction << std::endl;
    n_blocks = 8;
    block_ptr = 0;
    accepted = std::vector <double>(n_params);
    rejected = std::vector <double>(n_params);
    proposal_sd = input_var;
    llik = 0.0;
    lprior = 0.0;
    transitions = arma::zeros(n_tot+1, n_tot+1);
    stationary_prev = arma::zeros(n_tot+1, n_tot+1);
    transitions_t = arma::zeros(n_tot+1, n_tot+1);
    use_mean_ab = mean_ab;
    update_transitions();
}

void Param::update_transitions () {
    double tot_rate = 0.0;
    for (int i=1; i<=n_tot; i++) {
        transitions(i, 0) = mu[i-1];
        transitions(0, i) = lambda[i-1];
        tot_rate += transitions(0, i);
    }
    double curr_interaction;
    transitions(0, 0) = -tot_rate;
    for (int row_i=1; row_i<=n_tot; row_i++) {
        tot_rate = transitions(row_i, 0);
        for (int col_i=1; col_i<=n_tot; col_i++) {
            if (row_i!=col_i) {
                if (col_i <= n_vtypes) curr_interaction = interaction[0];
                else curr_interaction = interaction[1];
                transitions(row_i, col_i) = lambda[col_i-1]*curr_interaction;
                tot_rate += transitions(row_i, col_i);
            }
        }
        transitions(row_i, row_i) = -tot_rate;
    }
}

void Param::next_block() {
    block_ptr++;
//    while ((block_ptr >= n_tot) & (block_ptr < 2*n_tot)) block_ptr++;
    if (block_ptr >= n_params) block_ptr = 0;
}

void calc_expm(bool vacc, int ind, Data data, arma::mat& original_matrix, arma::mat& transitions1, double multiplier, int n_tot, int n_vtypes, std::vector<double> p0, std::vector <double> thetaSI, bool use_mean_ab) {
    int n_vacc = data["n_vacc"];
    arma::mat matrix_to_change(transitions1);
    double tot_rate = 0.0;
    double full_immunity;
    double susceptibility;
    if (use_mean_ab) susceptibility = data.get_mean_ab_ind(ind);
    for (int i=0; i<n_tot; i++) {
        matrix_to_change(0, i+1) = transitions1(0, i+1)*multiplier;
        matrix_to_change(i+1, 0) = transitions1(i+1, 0)*multiplier;
        if (vacc & (i<n_vtypes) & (multiplier < STATIONARY_TIME)) {
            full_immunity = 1.0;//(double) (ind >= (n_vacc*p0[i]));
            if (!use_mean_ab) susceptibility = data.get_ab(ind, i);
            //matrix_to_change(0, i+1) *= p0[i]-susceptibility*thetaSI[i];//*full_immunity*ind_frailty_SI[(i-1)*n_vacc+ind]
            matrix_to_change(0, i+1) *= p0[i]*std::exp(-susceptibility*thetaSI[i]);
            matrix_to_change(i+1, 0) *= 1.0;//thetaIS[i];//*full_immunity;//*ind_frailty_IS[ind];
        }
        if (matrix_to_change(0, i+1)<0.0) matrix_to_change(0, i+1) = 0.0;
        tot_rate += matrix_to_change(0, i+1);
    }
    matrix_to_change(0, 0) = -tot_rate;
    for (int row_i=1; row_i<=n_tot; row_i++) {
        tot_rate = matrix_to_change(row_i, 0);
        for (int col_i=1; col_i<=n_tot; col_i++) {
            if (row_i!=col_i) {
                matrix_to_change(row_i, col_i) = transitions1(row_i, col_i)*multiplier;
                if (vacc & (col_i<=n_vtypes) & (multiplier < STATIONARY_TIME)) {
                    full_immunity =1.0;// (double) (ind >= (n_vacc*p0[col_i-1]));
                    if (!use_mean_ab) susceptibility = data.get_ab(ind, col_i-1);
                    //matrix_to_change(row_i, col_i) *= p0[col_i-1]-susceptibility*thetaSI[col_i-1];//*full_immunity*ind_frailty_SI[(col_i-1)*n_vacc+ind]
                    matrix_to_change(row_i, col_i) *= p0[col_i-1]*std::exp(-susceptibility*thetaSI[col_i-1]);
                }
                if (matrix_to_change(row_i, col_i)<0.0) matrix_to_change(row_i, col_i) = 0.0;
                tot_rate += matrix_to_change(row_i, col_i);
            }
        }
        matrix_to_change(row_i, row_i) = -tot_rate;
    }
    try {
        original_matrix = arma::expmat(matrix_to_change);
//        print_matrix(original_matrix, "matrix_exp.txt", n_tot);
//        print_matrix(matrix_to_change, "matrix.txt", n_tot);
    }
    catch(...) {
        #pragma omp critical
        print_matrix(matrix_to_change, "matrix.txt", n_tot);
        std::cout << "EXPM failed" << std::endl;
    }
}

double Param::calc_llik (Data data) {
    int total_threads = omp_get_num_threads();
//    std::cout << total_threads << " threads" << std::endl;
    std::vector<double>newllik(total_threads, 0.0);
    std::vector <Data> data_vec;
    std::vector <arma::mat> transitions_t_vec;
    std::vector <arma::mat> stationary_prev_vec;
    std::vector <arma::mat> transitions_vec;
    std::vector <std::vector <double> > p0_vec;
    std::vector <std::vector <double> > thetaSI_vec;
    for (int i=0; i<total_threads; ++i) {
        data_vec.push_back(data);
        transitions_t_vec.push_back(transitions_t);
        transitions_vec.push_back(transitions);
        stationary_prev_vec.push_back(stationary_prev);
        p0_vec.push_back(p0);
        thetaSI_vec.push_back(thetaSI);
    }
#pragma omp parallel for schedule(static, 1)
    for (int tn=0; tn<total_threads; tn++) {
        for (int i=tn; i<data_vec[tn]["n_vacc"]; i+=total_threads) {
//        for (int i=tn; i<200; i+=total_threads) {
            int position1, position2;
//            tn = omp_get_thread_num();
            /*For each individual, calculate the likelihood of observations*/
            position1 = (int) data_vec[tn].get_swab_v(i, 0);
            // calculate the probability of first swab results, i.e. at stationarity
            try{
                calc_expm(true, i, data_vec[tn], stationary_prev_vec[tn], transitions_vec[tn], STATIONARY_TIME, n_tot, n_vtypes, p0_vec[tn], thetaSI_vec[tn], use_mean_ab);
                newllik[tn] += log(stationary_prev_vec[tn](0, position1));
            }
            catch (...) {
#pragma omp critical
                newllik[tn] += SMALLEST_NUMBER;
                break;
                //return (SMALLEST_NUMBER);
            }
            for (int time_step=1; time_step<data["n_swabs"]; time_step++) {
//                if (tn==0) std::cout << i << " " << time_step << " " << newllik[tn] << " prob: " << stationary_prev_vec[tn](0, position1) << " Ab: " << data_vec[tn].get_ab(i, 0)<< std::endl;
                try{
                    calc_expm(true, i, data_vec[tn], transitions_t_vec[tn], transitions_vec[tn], data_vec[tn].get_swab_time(time_step), n_tot, n_vtypes, p0_vec[tn], thetaSI_vec[tn], use_mean_ab);
                }
                catch (...) {
#pragma omp critical
                    newllik[tn] += SMALLEST_NUMBER;
                    break;
//                    return (SMALLEST_NUMBER);
                }
                position2 = (int) data_vec[tn].get_swab_v(i, time_step);
                newllik[tn] += log(transitions_t_vec[tn](position1,position2));
                position1 = position2;
            }
        }
    }
//    try{
//        calc_expm(false, 0, data, stationary_prev, STATIONARY_TIME);
//    }
//    catch (...) {
//        return (SMALLEST_NUMBER);
//    }
//    for (int i=0; i<data["n_nvacc"]; i++) {
//        position1 = (int) data.get_swab_nv(i, 0);
//        newllik += log(stationary_prev(0, position1));
//        for (int time_step=1; time_step<data["n_swabs"]; time_step++) {
//            try{
//                calc_expm(false, i, data, transitions_t, data.get_swab_time(time_step));
//            }
//            catch (...) {
//                return (SMALLEST_NUMBER);
//            }
//            position2 = (int) data.get_swab_nv(i, time_step);
//            newllik += log(transitions_t(position1, position2));
//            position1 = position2;
//        }
//    }
//    if (!std::isfinite(newllik)) newllik = SMALLEST_NUMBER;
    double sum_newllik = 0.0;
    for (int i=0; i!=newllik.size(); ++i) {
        if(!std::isfinite(newllik[i])) return (SMALLEST_NUMBER);
        else if (newllik[i]==SMALLEST_NUMBER) return (SMALLEST_NUMBER);
        else sum_newllik += newllik[i];
    }
    return (sum_newllik);
}

double Param::calc_lprior (int block_i) const {
    double newlprior = 0.0;
    if (block_i<n_tot) {
//        for (auto i=lambda.begin(); i<lambda.end(); i++) {
            boost::math::uniform_distribution <double> density(0.000000001, 1.0);
            newlprior += log(boost::math::pdf(density, lambda[block_i]));
//        }
    }
    else if (block_i<n_tot*2) {
//        for (auto i=mu.begin(); i<mu.end(); i++) {
            boost::math::uniform_distribution <double> density(0.000000001, 0.05);
            newlprior += log(boost::math::pdf(density, mu[block_i-n_tot]));
//        }
    }
    else if (block_i<n_tot*2+n_vtypes) {
//        for (auto i=thetaSI.begin(); i<thetaSI.end(); i++) {
            boost::math::uniform_distribution<>density(-1.0, 5.0);
            newlprior += log(boost::math::pdf(density, thetaSI[block_i-n_tot*2]));
//        }
    }
//    else if (block_i==3) {
//        for (auto i=thetaIS.begin(); i<thetaIS.end(); i++) {
//            boost::math::uniform_distribution<>density(0.5, 2.0);
//            newlprior += log(boost::math::pdf(density, *i));
//        }
//    }
    else if (block_i<n_tot*2+n_vtypes*2) {
//        for (auto i=p0.begin(); i<p0.end(); i++) {
            boost::math::uniform_distribution<>density(0.0, 2.0);
            newlprior += log(boost::math::pdf(density, p0[block_i-n_tot*2-n_vtypes]));
//        }
    }
//    else if (block_i==5) {
//        boost::math::exponential_distribution<double>density(1.0);
//        newlprior += log(boost::math::pdf(density, frailtySI));
//    }
//    else if (block_i==6) {
//        boost::math::exponential_distribution<double>density(1.0);
//        newlprior += log(boost::math::pdf(density, frailtyIS));
//    }
    else if (block_i<n_tot*2+n_vtypes*2+2) {
        boost::math::beta_distribution<double>density(200.0, 300.0);// mean = 0.4, sd = 0.2
        if (interaction[block_i-n_tot*2-2*n_vtypes] < 0.0 | interaction[block_i-n_tot*2-2*n_vtypes] > 1.0) newlprior = SMALLEST_NUMBER;
        else newlprior += log(boost::math::pdf(density, interaction[block_i-n_tot*2-2*n_vtypes]));
    }
    if (!std::isfinite(newlprior)) newlprior = SMALLEST_NUMBER;
    return (newlprior);
}

double Param::calc_lprior () const {
    double newlprior = 0.0;
    for (int i=0; i<n_params; i++) {
        newlprior += calc_lprior(i);
        if (newlprior <= SMALLEST_NUMBER) return (SMALLEST_NUMBER);
    }
    return (newlprior);
}

double Param::uni_propose (double oldval, double sd, int ntries) const {
    boost::normal_distribution<> nd(oldval, sd);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(rng, nd);
    double newval = var_nor();
    if (ntries > 0) {
        while (newval < 0) {
            newval = var_nor();
            ntries--;
            if (ntries < 0) break;
        }
    }
    return (newval);
}

void Param::alter_param (bool reject) {
    int ntries = 10;
    int ptr;
    if (block_ptr<n_tot) {
        if (reject) {
            lambda = tempparam;
            update_transitions();
            return;
        }
        tempparam = lambda;
        //for (int i=0; i<n_tot; i++) {
        lambda[block_ptr] = uni_propose(tempparam[block_ptr], proposal_sd[block_ptr], ntries);
        //}
        update_transitions();
    }
    else if (block_ptr<n_tot*2) {
        if (reject) {
            mu = tempparam;
            update_transitions();
            return;
        }
        tempparam = mu;
        //for (int i=0; i<n_tot; i++) {
        ptr = block_ptr-n_tot;
            mu[ptr] = uni_propose(tempparam[ptr], proposal_sd[block_ptr], ntries);
        //}
        update_transitions();
    }
    else if (block_ptr<n_tot*2+n_vtypes) {
        if (reject) {
            thetaSI = tempparam;
            return;
        }
        tempparam = thetaSI;
//        for (int i=0; i<n_vtypes; i++) {
            ptr = block_ptr - n_tot*2;
            thetaSI[ptr] = uni_propose(tempparam[ptr], proposal_sd[block_ptr], 0);
//        }
    }
//    else if (block_ptr==3) {
//        if (reject) {
//            thetaIS = tempparam;
//            return;
//        }
//        tempparam = thetaIS;
//        for (int i=0; i<n_vtypes; i++) {
//            thetaIS[i] = uni_propose(tempparam[i], proposal_sd[block_ptr], ntries);
//        }
//    }
    else if (block_ptr<n_tot*2+n_vtypes*2) {
        if (reject) {
            p0 = tempparam;
            return;
        }
        tempparam = p0;
//        for (int i=0; i<n_vtypes; i++) {
        ptr = block_ptr - (n_tot*2+n_vtypes);
            p0[ptr] = uni_propose(tempparam[ptr], proposal_sd[block_ptr], 0);
//        }
    }
//    else if (block_ptr==5) {
//        if (reject) {
//            frailtySI = tempparam[0];
//            return;
//        }
//        tempparam[0] = frailtySI;
//        frailtySI = uni_propose(tempparam[0], proposal_sd[block_ptr], ntries);
//    }
//    else if (block_ptr==6) {
//        if (reject) {
//            frailtyIS = tempparam[0];
//            return;
//        }
//        tempparam[0] = frailtyIS;
//        frailtyIS = uni_propose(tempparam[0], proposal_sd[block_ptr], ntries);
//    }
    else if (block_ptr<n_tot*2+n_vtypes*2+2) {
        ptr = block_ptr - n_tot*2 - n_vtypes*2;
        if (reject) {
            interaction[ptr] = tempparam[ptr];
            update_transitions();
            return;
        }
        tempparam[ptr] = interaction[ptr];
        interaction[ptr] = uni_propose(tempparam[ptr], proposal_sd[block_ptr], ntries);
        update_transitions();
    }
}

void Param::mcmc_move (Data data, bool adapt, double optimal_adapt) {
    next_block();
    // Code to propose
    double previous_val = (*this)[block_ptr];
    alter_param (false);
    get_rand_frailty(data);
    double newlprior = calc_lprior();
    double newllik = SMALLEST_NUMBER;
    if (newlprior > SMALLEST_NUMBER) newllik = calc_llik(data);
    boost::uniform_01<> zeroone;
    boost::variate_generator<boost::mt19937&, boost::uniform_01<> > runif(rng, zeroone);
    double z = log(runif());
    double ratio = (newllik+newlprior-llik-lprior);
//    boost::math::normal_distribution<double>norm(0.0,1.0);
//    double temp1 = boost::math::cdf(norm, -previous_val/proposal_sd[block_ptr]);
//    double temp2 = boost::math::cdf(norm, -(*this)[block_ptr]/proposal_sd[block_ptr]);
//    double K = (1.0-temp1)/(1.0-temp2);
//    ratio += log(K);
    bool reject = z > ratio;
    if (reject) {
        alter_param (true);
        rejected[block_ptr]++;
    } else {
        accepted[block_ptr]++;
        llik = newllik;
        lprior = newlprior;
    }
    if (adapt) {
        double acceptance_rate = (double)accepted[block_ptr] / (double)(accepted[block_ptr]+rejected[block_ptr]);
        double change = exp(0.999/2.0*(acceptance_rate-optimal_adapt));
        proposal_sd[block_ptr] *= change;
        std::fill(accepted.begin(), accepted.end(), 0.0);
        std::fill(rejected.begin(), rejected.end(), 0.0);
    }
}


double Param::operator[](int i) {
    double value;
    if (i < n_tot) value = lambda[i];
    else if (i < n_tot*2) value = mu[i-n_tot];
    else if (i < n_tot*2+n_vtypes) value = thetaSI[i-n_tot*2];
//    else if (i < n_tot*2+n_vtypes*2) value = thetaIS[i-n_tot*2-n_vtypes];
    else if (i < n_tot*2+n_vtypes*2) value = p0[i-n_tot*2-n_vtypes];
//    else if (i == n_tot*2+n_vtypes*3) value = frailtySI;
//    else if (i == n_tot*2+n_vtypes*3+1) value = frailtyIS;
    else value = interaction[i-n_tot*2-n_vtypes*2];
    return (value);
}

void Param::get_rand_frailty (Data & data) {
    ind_frailty_SI.clear();
    // ind_frailty_IS.clear();
    for (int type_i=0; type_i<data["n_vtypes"]; type_i++) {
        double m1 = data.get_mean_ab(type_i);
        for (int ind_i=0; ind_i<data["n_vacc"]; ind_i++) {
            // ind_frailty_SI.push_back((data.get_ab(ind_i, type_i)*sqrt(frailtySI)/sd1 + (1.0-m1)));
            ind_frailty_SI.push_back((data.get_max_ab(type_i)-data.get_ab(ind_i, type_i))/m1);
        }
        // ind_frailty_SI.push_back(R::rgamma(1.0/frailtySI, frailtySI));
        // ind_frailty_IS.push_back(R::rgamma(1.0/frailtyIS, frailtyIS));
    }
}

void Param::initial_calc(Data data) {
    // get_rand_frailty(data);
    lprior = calc_lprior();
    if (lprior > SMALLEST_NUMBER) llik = calc_llik(data);
    else llik = SMALLEST_NUMBER;
}

void Param::initialize_file (std::string filename) {
    std::ofstream output_file;
    output_file.open(filename);
    output_file << "state\tposterior\tlikelihood\tprior";
    for (int i=0; i<n_tot; i++) output_file << "\tlambda" << i;
    for (int i=0; i<n_tot; i++) output_file << "\tmu" << i;
    for (int i=0; i<n_vtypes; i++) output_file << "\tthetaSI" << i;
//    for (int i=0; i<n_vtypes; i++) output_file << "\tthetaIS" << i;
    for (int i=0; i<n_vtypes; i++) output_file << "\tp0" << i;
    for (int i=0; i<2; i++) output_file << "\tinteraction" << i;
//    output_file <<"\tfrailtySI\tfrailtyIS\tinteraction";
    output_file << std::endl;
    output_file.close();
}

void Param::print_to_file (std::string filename, int iter, arma::mat &output_mat) {
    output_mat(iter, 0) = iter;
    output_mat(iter, 1) = llik+lprior;
    output_mat(iter, 2) = llik;
    output_mat(iter, 3) = lprior;
    std::ofstream output_file;
    output_file.open(filename, std::ofstream::app);
    output_file << iter << "\t" << llik+lprior << "\t" << llik << "\t" << lprior;
    for (int i=0; i<n_params; i++) {
        output_file << "\t" << (*this)[i];
        output_mat(iter, i+4) = (*this)[i];
    }
    output_file << std::endl;
    output_file.close();
}

