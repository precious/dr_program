#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

typedef float real;
typedef float velocity;

extern velocity ORBITAL_VELOCITY;
//extern velocity ION_VELOCITY;
//extern velocity ELECTRON_VELOCITY;

extern real ELECTRON_VELOCITY_M;
extern real ELECTRON_VELOCITY_D;
extern real ION_VELOCITY_M;
extern real ION_VELOCITY_D;

// заряд, умноженный на 10^10
extern double ELECTRON_ELECTRIC_CHARGE;
extern double ION_ELECTRIC_CHARGE;

// Дебаевский радиус
extern int ELECTRONS_GENERATIVE_SPHERE_RADIUS;
extern int IONS_GENERATIVE_SPHERE_RADIUS;

// плотность частиц
extern int ELECTRONS_CONSISTENCE;
extern int IONS_CONSISTENCE;

#endif // CONSTANTS_H

