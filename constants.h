#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

#define EXIT_ERR(msg) { cerr << msg << "\nerrno: " << errno << endl; Graphics::quitGraphics(1); }
#define PRINTLN(arg) cout << arg << endl
#define PRINT(arg) cout << arg && cout.flush()
#define COUT(args) cout << args << endl

// Particle types
#define PTYPE_ELECTRON 1
#define PTYPE_ION 0

// charges
extern double ELECTRON_ELECTRIC_CHARGE;
extern double ION_ELECTRIC_CHARGE;

// charge to mass ratio
extern double ELECTRON_CHARGE_TO_MASS;
extern double ION_CHARGE_TO_MASS;

// current density
extern double ELECTRON_CURRENT_DENSITY;
extern double ION_CURRENT_DENSITY;

#define PARTICLE_CHARGE(_particle_type) \
    ((_particle_type == PTYPE_ELECTRON)? ELECTRON_ELECTRIC_CHARGE: ION_ELECTRIC_CHARGE)

#define PARTICLE_CHARGE_TO_MASS(_particle_type) \
    ((_particle_type == PTYPE_ELECTRON)? ELECTRON_CHARGE_TO_MASS: ION_CHARGE_TO_MASS)

#define PARTICLE_CURRENT_DENSITY(_particle_type) \
    ((_particle_type == PTYPE_ELECTRON)? ELECTRON_CURRENT_DENSITY: ION_CURRENT_DENSITY)

typedef double real;
typedef float velocity;

extern velocity ORBITAL_VELOCITY;
//extern velocity ION_VELOCITY;
//extern velocity ELECTRON_VELOCITY;

extern real ELECTRON_VELOCITY_M;
extern real ELECTRON_VELOCITY_D;
extern real ION_VELOCITY_M;
extern real ION_VELOCITY_D;

// Debye radius
extern int ELECTRONS_GENERATIVE_SPHERE_RADIUS;
extern int IONS_GENERATIVE_SPHERE_RADIUS;

// particles density
extern int ELECTRONS_CONSISTENCE;
extern int IONS_CONSISTENCE;

extern double VACUUM_PERMITTIVITY;
// e0 ~ 8.854187817620*10^-12 F/m

#endif // CONSTANTS_H

