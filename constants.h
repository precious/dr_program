#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

#define EXIT_ERR(msg) { cerr << msg << "\nerrno: " << errno << endl; Graphics::quitGraphics(1); }
#define PRINTLN(arg) cout << arg << endl
#define PRINT(arg) cout << arg && cout.flush()
#define COUT(args) cout << args << endl

typedef double real;
typedef float velocity;

extern velocity ORBITAL_VELOCITY;
//extern velocity ION_VELOCITY;
//extern velocity ELECTRON_VELOCITY;

extern real ELECTRON_VELOCITY_M;
extern real ELECTRON_VELOCITY_D;
extern real ION_VELOCITY_M;
extern real ION_VELOCITY_D;

// charge multiplied by 10^10
extern double ELECTRON_ELECTRIC_CHARGE;
extern double ION_ELECTRIC_CHARGE;

// Debye radius
extern int ELECTRONS_GENERATIVE_SPHERE_RADIUS;
extern int IONS_GENERATIVE_SPHERE_RADIUS;

// particles density
extern int ELECTRONS_CONSISTENCE;
extern int IONS_CONSISTENCE;

// current density
extern double ELECTRON_CURRENT_DENSITY;
extern double ION_CURRENT_DENSITY;

extern double VACUUM_PERMITTIVITY;
// e0 ~ 8.854187817620*10^-12 F/m

#endif // CONSTANTS_H

