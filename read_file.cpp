#include "read_file.h"

vector<PlaneType>* getCoordinatesFromFile(char *filename) {
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
            } else {
                tempThreePoints->set[i].x *= SCALE;
                tempThreePoints->set[i].y *= SCALE;
                tempThreePoints->set[i].z *= SCALE;
            }
        }

        // if all values have been read successfuly then push array to result list
        if (i == 3)
            coordinatesList->push_back(PlaneType(*tempThreePoints));
        // then delete it
        delete tempThreePoints;
        // scipping remaining characters in current string
        while (!fileInputStream.eof() && fileInputStream.get() != '\n');
    }

    fb.close();
    return coordinatesList;
}
