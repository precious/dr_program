#include "constants.h"

velocity ORBITAL_VELOCITY = 7907.343098064;

real ELECTRONS_GENERATIVE_SPHERE_RADIUS = 131.414769518;
real IONS_GENERATIVE_SPHERE_RADIUS = ELECTRONS_GENERATIVE_SPHERE_RADIUS;

real ELECTRON_VELOCITY_M = 110901445.521404395; // 2*sqrt(2*1.6*10^-19*27.5*10^3/(9.11*10^-31*pi))
real ELECTRON_VELOCITY_D = 69497174.760350449; // sqrt(1.6*10^-19*27.5*10^3/(9.11*10^-31))
real ION_VELOCITY_M = 2613670.454744482; // 2*sqrt(2*1.6*10^-19*28*10^3/(1.67*10^-27*pi))
real ION_VELOCITY_D = 1637875.065607546; // sqrt(1.6*10^-19*28*10^3/(1.67*10^-27))


double ELECTRON_ELECTRIC_CHARGE = -0.000000000000000000160217656535;
double ION_ELECTRIC_CHARGE = 0.000000000000000000160217656535;

int ELECTRONS_CONSISTENCE = 1200000; // Novik, page 20
int IONS_CONSISTENCE = ELECTRONS_CONSISTENCE/2;//1300000; // Novik, page 20

double ELECTRON_CURRENT_DENSITY = -0.000005334;
// 1.60217656535*10^-19*1.2*10^6*sqrt(27.5*10^3*1.60217648740*10^-19/(2*pi*9.1093829140*10^-31));
double ION_CURRENT_DENSITY = 0.000000136;
// 1.60217656535*10^-19*1.3*10^6*sqrt(28*10^3*1.60217648740*10^-19/(2*pi*1.67262177774*10^-27))

double ELECTRON_CHARGE_TO_MASS = -175882008745.910974667;
// (-1.60217656535*10^-19)/(9.1093829140*10^-31)
double ION_CHARGE_TO_MASS = 95788335.813420796;
// (1.60217656535*10^-19)/(1.67262177774*10^-27)

double VACUUM_PERMITTIVITY = 0.00000000000885418782;

// E ~ q, proofs:
// http://ens.tpu.ru/POSOBIE_FIS_KUSN/%DD%EB%E5%EA%F2%F0%EE%F1%F2%E0%F2%E8%EA%E0.%20%CF%EE%F1%F2%EE%FF%ED%ED%FB%E9%20%D2%EE%EA/05-2.htm
// http://physic.kemsu.ru/pub/library/learn_pos/UMK_Electrostat/Pages/%D0%A2%D0%B5%D0%BE%D1%80%D0%B8%D1%8F/%D0%9F%D1%80%D0%BE%D0%B2%D0%BE%D0%B4%D0%BD%D0%B8%D0%BA%D0%B8%20%D0%B2%20%D1%8D%D0%BB%D0%B5%D0%BA%D1%80%D0%B8%D1%87%D0%B5%D1%81%D0%BA%D0%BE%D0%BC%20%D0%BF%D0%BE%D0%BB%D0%B5.htm
