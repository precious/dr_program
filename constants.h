#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <map>
#include <vector>
#include <string>

#define ORBITAL_UNIT 1
#define SI_UNIT 2
#define GL_UNIT 3

using namespace std;

typedef float velocity;

extern map<string,velocity> constantsMap;
velocity getVelocity(char*,int);
velocity getVelocityOrb(char*);
velocity getVelocitySi(char*);
velocity getVelocityGl(char*);
velocity convertVelocity(velocity,int,int);


#endif // CONSTANTS_H

