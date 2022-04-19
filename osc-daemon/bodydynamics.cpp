#include <iostream>
#include <vector>
#include <math.h>

#include "osctypes.hpp"
#include "planet.cpp"

namespace osc{

    rotStates rotationDynamicModel(rotStates curStates, rotStates control, std::array<std::array<double, 3>, 3> inverseInertia) { //someone make matrices inversible pls
    rotStates dotStates;
        dotStates.omega.data[0] = inverseInertia[0][0]*control.omega.data[0]+inverseInertia[0][1]*control.omega.data[1]+inverseInertia[0][2]*control.omega.data[2];
        dotStates.omega.data[1] = inverseInertia[1][0]*control.omega.data[0]+inverseInertia[1][1]*control.omega.data[1]+inverseInertia[1][2]*control.omega.data[2];
        dotStates.omega.data[2] = inverseInertia[2][0]*control.omega.data[0]+inverseInertia[2][1]*control.omega.data[1]+inverseInertia[2][2]*control.omega.data[2];
        dotStates.q             = curStates.q.operator/(2);
    };

    posStates positiondynamicmodel(posStates curStates, double thrust, double Ve) { //work in progress
    double r = curStates.r.mag();
    double v = curStates.v.mag();
    double C =  (3 * planet.J2 * planet.sgp * pow2(planet.sMa)) / (2 * pow(r, 5)); //constant value
    double accRelConst = (-planet.sgp / pow(r, 1.5));

        vec3 accRel = curStates.r.operator*(accRelConst);

         vec3 accGrav;
            accGrav.data[0] = C * (5 * pow2(curStates.r.data[2] / r) - 1) * curStates.r.data[0];
            accGrav.data[1] = C * (5 * pow2(curStates.r.data[2] / r) - 1) * curStates.r.data[1];
            accGrav.data[2] = C * (5 * pow2(curStates.r.data[2] / r) - 3) * curStates.r.data[2];

        vec3 accThrust = curStates.v.operator*(thrust / curStates.m / v);

        posStates dotStates;
            dotStates.r = curStates.v;
            dotStates.v = accRel.operatorplus2(accGrav, accThrust);
            dotStates.m = -thrust / Ve * curStates.m;
    };

}