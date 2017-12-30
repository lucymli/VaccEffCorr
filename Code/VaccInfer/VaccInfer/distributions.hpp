//
//  distributions.hpp
//  VaccInfer
//
//  Created by Lucy Li on 5/17/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#ifndef distributions_hpp
#define distributions_hpp

#include <iostream>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions.hpp>

double rnorm(double, double, double, double, int);

double rnorm(double, double);

double dnorm(double, double, double);

double runif(double, double);

double runif();

double dunif(double, double, double);

int rsample(int, int);

int rsample(std::vector <double>);

int rsample(int, std::vector <double>);

double rsample(std::vector <double>, std::vector <double>);

double get_rand_num(std::string, double, double);

double get_density(double, std::string, double, double, bool);

#endif /* distributions_hpp */
