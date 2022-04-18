#include <vector>
#include <array>
#include <math.h>

//default values for Earth
struct celestial{
    double sgp; //standard gravitational parameter mu
    double sMa; //equatorial radius of planet
    double flat; //ellipticity of planet's surface
    double sma; //polar radius of planet
    double ecc; //ecccentricity of planet's surface
    double J2; //second zonal gravitational harmonic
    celestial(double initsgp, double initsMa, double initfla, double initsma, double initecc, double initJ2):
        sgp(initsgp), sMa(initsMa), flat(initfla), sma(initsma), ecc(initecc), J2(initJ2){};
};
celestial earth(3.986004418e14, 6378137, 1/298.257223563,   6356752.314,    0.08181919084,  1082.63e-6);
celestial venus(3.24859e14,     6051800, 0,                 6051800,        0,              4.458e-6);
celestial mars(4.282837e13,     3396200, 0.00589,           3376200,        0.1083757717,   1960.45e-6);
celestial moon(4.9048695e12,    1738100, 0.0012,            1736000,        0.04897509571,  202.7e-6);

// add more planets if we have time or if i can be bothered
celestial planet=earth;//initialisation value, have this in start-up of final code
//might need copy constructor here


