#include "time_utils.h"

int CLOCK_ID = CLOCK_THREAD_CPUTIME_ID;

void Time::printTimespec(timespec *ts) {
    cout << "seconds: " << ts->tv_sec << ", nanoseconds: " << ts->tv_nsec << endl;
}

timespec* Time::getTimespecDelta(timespec *older,timespec *newer) {
    newer->tv_nsec -= older->tv_nsec;
    newer->tv_sec -= older->tv_sec;
    if (newer->tv_nsec < 0) {
        newer->tv_sec--;
        newer->tv_nsec += 1000000000;
    }
    return newer;
}

double Time::getRandom() {
    static random_device seed;
    static UniformDistributionGenerator* generator =
            new UniformDistributionGenerator(Engine(seed()),UniformDistribution(0,1));
    return (*generator)();
}

double Time::getRandom(double left,double right) {
    return getRandom()*(right - left) + left;
}

UniformDistributionGenerator Time::getUniformDistributionGenerator(double min,double max) {
    static random_device seed;
    return UniformDistributionGenerator(Engine(seed()),UniformDistribution(min,max));
}

GaussianDistributionGenerator Time::getGaussianDistributionGenerator(double M,double D) {
    static random_device seed;
    return GaussianDistributionGenerator(Engine(seed()),GaussianDistribution(M,D));
}

MaxwellDistributionSpeedGenerator Time::getMaxwellDistributionSpeedGenerator(double M,double D) {
    return  [M,D]() -> velocity {
            /// TODO check for M, D
            static GaussianDistributionGenerator gdg = getGaussianDistributionGenerator(M/sqrt(3.0),D/sqrt(3.0));
            return sqrt(pow(gdg(),2) + pow(gdg(),2) + pow(gdg(),2));
    };
}

