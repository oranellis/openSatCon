#include <iostream>
#include <vector>
#include <math.h>

#include "osctypes.hpp"
#include "planet.cpp"

namespace osc{

    rotstates rotationdynamicmodel(rotstates curstates, rotstates control, std::array<std::array<double, 3>, 3> inverseinertia) { //someone make matrices inversible pls
    rotstates dotstates;
        dotstates.o1=inverseinertia[0][0]*control.o1+inverseinertia[0][1]*control.o2+inverseinertia[0][2]*control.o3;
        dotstates.o2=inverseinertia[1][0]*control.o1+inverseinertia[1][1]*control.o2+inverseinertia[1][2]*control.o3;
        dotstates.o3=inverseinertia[2][0]*control.o1+inverseinertia[2][1]*control.o2+inverseinertia[2][2]*control.o3;
        dotstates.q1=0.5*curstates.q1;
        dotstates.q2=0.5*curstates.q2;
        dotstates.q3=0.5*curstates.q3;
    };

    posstates positiondynamicmodel(posstates curstates, double thrust, double Ve) { //work in progress
    double r = sqrt(pow2(curstates.i)+pow2(curstates.j)+pow2(curstates.k));
    double v = sqrt(pow2(curstates.vi)+pow2(curstates.vj)+pow2(curstates.vk));
    double C =  (-3*planet.J2*planet.sgp*pow2(planet.sMa))/(2*pow(r,5))
                    *(1-(5*pow2(curstates.k)/pow2(r))); //constant value
        eci accgrav;
            accgrav.i=C*curstates.i;
            accgrav.j=C*curstates.j;
            accgrav.k=C*curstates.k;
        posstates dotstates;
            dotstates.i=(curstates.vi)*curstates.i;
            dotstates.j=(curstates.vj)*curstates.j;
            dotstates.k=(curstates.vk)*curstates.k;
            dotstates.vi=(-planet.sgp*curstates.i/pow3(r)+accgrav.i+((thrust*curstates.vi)/(curstates.m*v)))*curstates.vi;
            dotstates.vj=(-planet.sgp*curstates.j/pow3(r)+accgrav.j+((thrust*curstates.vj)/(curstates.m*v)))*curstates.vj;
            dotstates.vk=(-planet.sgp*curstates.k/pow3(r)+accgrav.k+((thrust*curstates.vk)/(curstates.m*v)))*curstates.vk;
            dotstates.m=-thrust/Ve*curstates.m;
    };

}