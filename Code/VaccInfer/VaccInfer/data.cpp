//
//  data.cpp
//  VaccInfer
//
//  Created by Lucy Li on 2/13/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#include <iostream>
#include "data.hpp"

Data::Data () {
    n_ind = 0;
    n_tot = 0;
    n_time = 0;
    n_vtypes = 0;
    n_predictors = 0;
}

Data::Data (std::string filename) {
    std::ifstream input;
    input.open(filename);
    std::string name;
    double number;
    int integer;
    while (!input.eof()) {
        input >> name;
        if (name == "n_ind") input >> n_ind;
        if (name == "n_tot") input >> n_tot;
        if (name == "n_time") input >> n_time;
        if (name == "n_predictors") input >> n_predictors;
        n_vtypes = n_predictors;
        if (name == "times") {
            for (int i=0; i<n_time; i++) {
                input >> number;
                times.push_back(number);
            }
        }
        if (name == "carriage") {
            carriage.resize(n_ind*n_time);
            for (int i=0; i<n_ind; i++) {
                for (int j=0; j<n_time; j++) {
                    input >> integer;
                    set_carriage(i, j, integer);
                }
            }
        }
        if (name == "predictors") {
            predictors.resize(n_ind*n_predictors);
            for (int i=0; i<n_ind; i++) {
                for (int j=0; j<n_predictors; j++) {
                    input >> number;
                    set_predictor(i, j, number);
                }
            }
        }
        if (name == "predictor_map") {
            predictor_map.resize(n_ind*n_predictors);
            for (int i=0; i<n_ind; i++) {
                for (int j=0; j<n_predictors; j++) {
                    input >> integer;
                    set_predictor_index(i, j, integer);
                }
            }
            break;
        }
    }
    input.close();
}

Data::Data (int num_ind, int num_types, int num_times, std::vector <double> time_points, int num_predictors) {
    // Generates an empty Data class
    for (int i=0; i<(num_ind*num_times); i++) carriage.push_back(0);
    for (int i=0; i<(num_ind*num_predictors); i++) {
        predictors.push_back(0);
        predictor_map.push_back(i%num_predictors);
    }
    for (int i=0; i<num_times; i++) {
        times.push_back(time_points[i]);
    }
    n_ind = num_ind;
    n_tot = num_types;
    n_time = num_times;
    n_vtypes = num_predictors;
    n_predictors = num_predictors;
}

int Data::get_carriage(int ind_i, int time_i) const {
    return (carriage[time_i * n_ind + ind_i]);
}

void Data::set_carriage(int ind_i, int time_i, int val) {
    carriage[time_i*n_ind + ind_i] = val;
}

double Data::get_metadata(int ind_i, int i) const {
    return (metadata[i*n_ind + ind_i]);
}

double Data::get_predictor(int ind_i, int predict_i) const {
    return (predictors[predict_i * n_ind + ind_i]);
}

void Data::set_predictor(int ind_i, int predictor_i, double val) {
    predictors[predictor_i*n_ind + ind_i] = val;
}

int Data::get_predictor_index(int ind_i, int predict_i) const {
    return (predictor_map[predict_i * n_ind + ind_i]);
}

void Data::set_predictor_index(int ind_i, int predictor_i, int mapping) {
    predictor_map[predictor_i*n_ind + ind_i] = mapping;
}

void Data::calc_mean_predictors() {
    for (int i=0; i<n_ind; i++) {
        mean_predictors.push_back(std::accumulate(predictors.begin()+i*n_predictors, predictors.begin()+(i+1)*n_predictors, 0.0));
        mean_predictors[n_ind] /= (double) n_predictors;
    }
}

void Data::print_data_to_file(std::string filename) {
    std::ofstream outputfile;
    outputfile.open(filename);
    outputfile << "n_ind " << n_ind << std::endl;
    outputfile << "n_tot " << n_tot << std::endl;
    outputfile << "n_time " << n_time << std::endl;
    outputfile << "n_predictors " << n_predictors << std::endl;
    outputfile << "times";
    for (int i=0; i<n_time; i++) outputfile << " " << times[i];
    outputfile << std::endl << "carriage" << std::endl;
    for (int i=0; i<n_ind; i++) {
        outputfile << get_carriage(i, 0);
        for (int j=1; j<n_time; j++) {
            outputfile << " " << get_carriage(i, j);
        }
        outputfile << std::endl;
    }
    outputfile << "predictors" << std::endl;
    for (int i=0; i<n_ind; i++) {
        outputfile << get_predictor(i, 0);
        for (int j=1; j<n_predictors; j++) {
            outputfile << " " << get_predictor(i, j);
        }
        outputfile << std::endl;
    }
    outputfile << "predictor_map" << std::endl;
    for (int i=0; i<n_ind; i++) {
        outputfile << get_predictor_index(i, 0);
        for (int j=1; j<n_predictors; j++) {
            outputfile << " " << get_predictor_index(i, j);
        }
        outputfile << std::endl;
    }
    outputfile.close();
}
