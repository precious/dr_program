#include "time_utils.h"

int CLOCK_ID = CLOCK_THREAD_CPUTIME_ID;

void printTimespec(timespec *ts) {
    cout << "seconds: " << ts->tv_sec << ", nanoseconds: " << ts->tv_nsec << endl;
}

timespec* getTimespecDelta(timespec *older,timespec *newer) {
    newer->tv_nsec -= older->tv_nsec;
    newer->tv_sec -= older->tv_sec;
    if (newer->tv_nsec < 0) {
        newer->tv_sec--;
        newer->tv_nsec += 1000000000;
    }
    return newer;
}

double getRandom() {
    static random_device seed;
    static UniformDistributionGenerator* generator =
            new UniformDistributionGenerator(Engine(seed()),UniformDistribution(0,1));
    return (*generator)();
}

UniformDistributionGenerator* getUniformDistributionGenerator(real min,real max) {
    static random_device seed;
    return new UniformDistributionGenerator(Engine(seed()),UniformDistribution(min,max));
}

GaussianDistributionGenerator* getGaussianDistributionGenerator(real M,real D) {
    static random_device seed;
    return new GaussianDistributionGenerator(Engine(seed()),GaussianDistribution(M,D));
}

