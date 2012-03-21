#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <map>
#include <vector>
#include <string>

using namespace std;

typedef float velocity;

extern map<string,velocity> constantsMap;
velocity getVelocity(char*);

#endif // CONSTANTS_H

