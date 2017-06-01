//
//  distributions.cpp
//  VaccInfer
//
//  Created by Lucy Li on 5/17/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#include "distributions.hpp"

boost::mt19937 rng(NULL);

double rnorm(double mean, double sd, double lower, double upper, int num_tries) {
    double val = rnorm(mean, sd);
    if (num_tries > 0) {
        while (val < lower | val > upper) {
            val = rnorm(mean, sd);
            num_tries--;
            if (num_tries < 0) break;
        }
    }
    return (val);
}

double rnorm(double mean, double sd) {
    boost::normal_distribution<> normal(mean, sd);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(rng, normal);
    double val = var_nor();
    return (val);
}

double dnorm(double x, double mean, double sd) {
    boost::math::normal_distribution<>density(mean, sd);
    double dens = log(boost::math::pdf(density, x));
    return (dens);
}

double runif(double lower, double upper) {
    boost::random::uniform_real_distribution <double> unif(lower, upper);
    double val = unif(rng);
    return (val);
}

double runif() {
    boost::uniform_01<> zeroone;
    boost::variate_generator<boost::mt19937&, boost::uniform_01<> > unif(rng, zeroone);
    double a = unif();
    return (a);
}

double dunif(double x, double lower, double upper) {
    boost::math::uniform_distribution <double> density (lower, upper);
    double dens = log(boost::math::pdf(density, x));
    return (dens);
}
