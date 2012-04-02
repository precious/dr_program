#ifndef READ_FILE_H
#define READ_FILE_H

#include <iostream>
#include <fstream>
#include <vector>
#include "types.h"
using namespace std;

typedef OrientedPlane PlaneType;

const float SCALE = 0.001; // coordinates are given in millimeters

vector<PlaneType>* getCoordinatesFromFile(char*);

#endif // READ_FILE_H
