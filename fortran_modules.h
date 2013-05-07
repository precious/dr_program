#ifndef FORTRAN_MODULES_H
#define FORTRAN_MODULES_H

#include "time_utils.h"

// the next funxtions are defined in Kul.f90 and linked from fortran_modules/Kul.o
extern "C" int laplace_(int*,Point*,Point*,Point*); // integer function Laplace(N,P1,P2,P3)
extern "C" void resultf_(Point*,real*,Vector*); //ResultF(Point, Pot, Grad)

int solveBoundaryProblem(vector<PlaneType> *coordinatesList,bool verbose = false) {
    // prepare arrays of vertices for algorithm
    verbose && COUT("initializing coordinates in fortran module...");
    int numVertices = coordinatesList->size();
    Point *P1 = new Point[numVertices];
    Point *P2 = new Point[numVertices];
    Point *P3 = new Point[numVertices];
    int i = 0;
    for(auto it = coordinatesList->begin();it != coordinatesList->end();++it,++i) {
        P1[i] = (*it).a;
        P2[i] = (*it).b;
        P3[i] = (*it).c;
    }

    // Solve boundary-valued problem
    timespec start, stop, *delta;
    verbose && COUT("solving boundary problem in fortran module...");
    clock_gettime(CLOCK_ID,&start);
    int retval = laplace_(&numVertices,P1,P2,P3);
    clock_gettime(CLOCK_ID,&stop);
    delta = Time::getTimespecDelta(&start,&stop);
    if (verbose) {
        PRINT("elapsed time: ");
        Time::printTimespec(delta);
    }
    delete[] P1;
    delete[] P2;
    delete[] P3;

    verbose && COUT("fortran init retval: " << retval);
    assert(retval > 0);
    return retval;
}

#endif // FORTRAN_MODULES_H
