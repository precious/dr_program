#ifndef READ_FILE_H
#define READ_FILE_H

#include <assimp/scene.h>
#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>

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
