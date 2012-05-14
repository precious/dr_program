#ifndef READ_FILE_H
#define READ_FILE_H

#include <assimp/assimp.hpp>
#include <assimp/aiScene.h>
#include <assimp/aiPostProcess.h>

#include <iostream>
#include <fstream>
#include <vector>
#include "types.h"

using namespace std;

typedef OrientedPlane PlaneType;

namespace File {
    extern float scaleFactor;

    vector<PlaneType>* getCoordinatesFromPlainFile(char*);
    vector<PlaneType>* getCoordinatesFromSpecialFile(char*);
}

#endif // READ_FILE_H
