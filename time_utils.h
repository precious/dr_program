#ifndef TIME_UTILS_H
#define TIME_UTILS_H

#include <unistd.h>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <random>
#include <functional>
#include "constants.h"

using namespace std;

#define RAND(n) rand()%n

extern int CLOCK_ID;

typedef uniform_real_distribution<double> UniformDistribution;
typedef normal_distribution<double> GaussianDistribution;

typedef mt19937 Engine;


template <typename Eng,typename Distrib>
class Generator {
private:
    Eng engine;
    Distrib distribution;
public:
    Generator(Eng e,Distrib d): engine(e), distribution(d) {}
    double operator() () {
        return distribution(engine);
    }
};

typedef Generator<Engine,UniformDistribution> UniformDistributionGenerator;
typedef Generator<Engine,GaussianDistribution> GaussianDistributionGenerator;
typedef function<velocity ()> MaxwellDistributionSpeedGenerator;

void printTimespec(timespec*);

UniformDistributionGenerator* getUniformDistributionGenerator(double,double);
GaussianDistributionGenerator* getGaussianDistributionGenerator(double,double);
MaxwellDistributionSpeedGenerator getMaxwellDistributionSpeedGenerator(double,double);

timespec* getTimespecDelta(timespec*,timespec*);

double getRandom();

double getRandom(double,double);

#endif // TIME_UTILS_H

