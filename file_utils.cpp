#include "file_utils.h"

float File::scaleFactor = 1; // by default coordinates are given in meters

vector<PlaneType>* File::getCoordinatesFromPlainFile(char *filename) {
    filebuf fb;
    if (!fb.open(filename,ios::in)) {
        cerr << "An error occurred while opening file" << endl;
        return NULL;
    }
    istream fileInputStream(&fb);
    vector<PlaneType>* coordinatesList = new vector<PlaneType>();
    ThreePoints *tempThreePoints;

    int i;
    while(!fileInputStream.eof()) {
        tempThreePoints = new ThreePoints();

        // read 3 points
        for (i = 0;i < 3;i++) {
            fileInputStream >> tempThreePoints->set[i].x
                            >> tempThreePoints->set[i].y
                            >> tempThreePoints->set[i].z;
            if (fileInputStream.fail()) {
                if (!fileInputStream.eof())
                    fileInputStream.clear();
                break;
            } /*else {
                tempThreePoints->set[i].x *= SCALE;
                tempThreePoints->set[i].y *= SCALE;
                tempThreePoints->set[i].z *= SCALE;
            }*/
        }

        // if all values have been read successfuly then push array to result list
        if (i == 3) {
            try {
                PlaneType *pt = new PlaneType(*tempThreePoints);
                coordinatesList->push_back(*pt);
                for (i = 0;i < 3;i++) {
                    coordinatesList->back().set[i].x *= scaleFactor;
                    coordinatesList->back().set[i].y *= scaleFactor;
                    coordinatesList->back().set[i].z *= scaleFactor;
                }
            } catch (ZeroNormal zn) {}

        }
        // then delete it
        delete tempThreePoints;
        // scipping remaining characters in current string
        while (!fileInputStream.eof() && fileInputStream.get() != '\n');
    }

    fb.close();
    return coordinatesList;
}

vector<PlaneType>* File::getCoordinatesFromSpecialFile(char *filename) {
    Assimp::Importer importer;
    aiScene *aiscene = (aiScene*)importer
        .ReadFile(filename,aiProcess_Triangulate|aiProcess_FixInfacingNormals|aiProcess_FindDegenerates
                  |aiProcess_PreTransformVertices|aiProcess_OptimizeMeshes|aiProcess_FindInvalidData|aiProcess_RemoveRedundantMaterials);
    if (aiscene == NULL) {
        cerr << "An error occurred while opening file" << endl;
        return NULL;
    }
    vector<PlaneType>* coordinatesList = new vector<PlaneType>();
    ThreePoints *tempThreePoints;
    int failedNumber = 0;

    for(unsigned int i = 0;i < aiscene->mNumMeshes;++i) {
        if (aiscene->mMeshes[i]->HasFaces()) {
            for(unsigned int j = 0;j < aiscene->mMeshes[i]->mNumFaces;++j) {

                tempThreePoints = new ThreePoints();

                for(int k = 0;k < 3;k++) {
                    tempThreePoints->set[k].x = aiscene->mMeshes[i]->mVertices[aiscene->mMeshes[i]->mFaces[j].mIndices[k]].x;
                    tempThreePoints->set[k].y = aiscene->mMeshes[i]->mVertices[aiscene->mMeshes[i]->mFaces[j].mIndices[k]].y;
                    tempThreePoints->set[k].z = aiscene->mMeshes[i]->mVertices[aiscene->mMeshes[i]->mFaces[j].mIndices[k]].z;
                }

                try {
                    PlaneType *pt = new PlaneType(*tempThreePoints);
                    coordinatesList->push_back(*pt);
                    for (int g = 0;g < 3;g++) {
                        coordinatesList->back().set[g].x *= scaleFactor;
                        coordinatesList->back().set[g].y *= scaleFactor;
                        coordinatesList->back().set[g].z *= scaleFactor;
                    }
                } catch (ZeroNormal zn) {
                    ++failedNumber;
                }

                delete tempThreePoints;
            }
        }
    }

    failedNumber && cerr << "failed: " << failedNumber << " polygons" << endl;
    return coordinatesList;
}
