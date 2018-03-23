//
//  distributions.cpp
//  VaccInfer
//
//  Created by Lucy Li on 5/17/17.
//  Copyright © 2017 Lucy M Li, CCDD, HSPH. All rights reserved.
//

#include "distributions.hpp"

boost::mt19937 rng(time(NULL));

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

double rlnorm(double mean, double sd, double lower, double upper, int num_tries) {
    double val = rlnorm(mean, sd);
    if (num_tries > 0) {
        while (val < lower | val > upper) {
            val = rlnorm(mean, sd);
            num_tries--;
            if (num_tries < 0) break;
        }
    }
    return (val);
}


double rlnorm(double mean, double sd) {
    boost::lognormal_distribution<> lognormal(mean, sd);
    boost::variate_generator<boost::mt19937&, boost::lognormal_distribution<> > var_lnor(rng, lognormal);
    double val = var_lnor();
    return (val);
}

double dlnorm(double x, double mean, double sd) {
    boost::math::lognormal_distribution<>density(mean, sd);
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

int rsample (int lower, int upper) {
    boost::random::uniform_int_distribution<> dist(lower, upper);
    return dist(rng);
}

int rsample(std::vector <double> probs) {
    boost::random::discrete_distribution<> dist(probs);
    return dist(rng);
}

int rsample(int lower, std::vector <double> probs) {
    int position = rsample(probs);
    return (lower + position);
}

double rsample(std::vector <double> numbers, std::vector <double> probs) {
    int position = rsample(probs);
    return (numbers[position]);
}

double get_density(double value, std::string distribution, double par1, double par2, bool return_log=true) {
    double density = 0.0;
    if (distribution == "unif") {
        density += dunif(value, par1, par2);
    }
    else if (distribution == "norm") {
        density += dnorm(value, par1, par2);
    }
    else if (distribution == "lnorm") {
        density += dlnorm(value, par1, par2);
    }
    if (!return_log) density = std::exp(density);
    return (density);
}

double get_rand_num(std::string distribution, double par1, double par2) {
    double rnum;
    if (distribution == "unif") {
        rnum = runif(par1, par2);
    }
    else if (distribution == "norm") {
        rnum = rnorm(par1, par2);
    }
    else if (distribution == "lnorm") {
        rnum = rlnorm(par1, par2);
    }
    return (rnum);
}
