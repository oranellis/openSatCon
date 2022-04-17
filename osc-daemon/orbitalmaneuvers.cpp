#include <iostream>
#include <vector>
#include <math.h>

#include "osctypes.hpp"
#include "planet.cpp"
#include "axistransforms.cpp"

namespace osc{

    double massburned(double dV, double mo, double Isp){
        double propuse;

        double mf = mo * exp(-dV/Isp);
        propuse=mo-mf;
        return propuse;
    };//calculate propellant mass used for a given delta V

};