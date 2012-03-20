#include "constants.h"

map<char*,velocity> *initVelocityMap() {
    map<char*,velocity> *vm = new map<char*,velocity>();
    vm->at("orbitalVelocity") = 7907.343098064;
    vm->at("ionVelocity") = 9090.364708796;
    vm->at("electronVelocity") = 389206.742021674;
    vm->at("maxElectronVelocity") = vm->at("electronVelocity");
    vm->at("minElectronVelocity") = vm->at("electronVelocity")*0.5;
    vm->at("maxIonVelocity") = vm->at("ionVelocity");
    vm->at("minIonVelocity") = vm->at("ionVelocity")*0.5;
}

map<char*,velocity> *velocityMap = initVelocityMap();

velocity getVelocitySi(char* velocityName) {
    return velocityMap->at(velocityName);
}

velocity getVelocityGl(char* velocityName) {
    return velocityMap->at(velocityName)*1000;
}
