#include <iostream>
#include <vector>
#include <math.h>

#include "osctypes.hpp"
#include "planet.cpp"

namespace osc{

    rotstates rotationdynamicmodel(rotstates curstates, rotstates control, inertiamatrix inverseinertia){//someone make matrices inversible pls
    rotstates dotstates;
        dotstates.o1=inverseinertia.J11*control.o1+inverseinertia.J12*control.o2+inverseinertia.J13*control.o3;
        dotstates.o2=inverseinertia.J21*control.o1+inverseinertia.J22*control.o2+inverseinertia.J23*control.o3;
        dotstates.o3=inverseinertia.J31*control.o1+inverseinertia.J32*control.o2+inverseinertia.J33*control.o3;
        dotstates.q1=0.5*curstates.q1;
        dotstates.q2=0.5*curstates.q2;
        dotstates.q3=0.5*curstates.q3;
    };

    posstates positiondynamicmodel(posstates curstates, posstates control){//work in progress
    double r = sqrt(pow(curstates.i,2)+pow(curstates.j,2)+pow(curstates.k,2));
    double C =  (-3*planet.J2*planet.sgp*pow(planet.sMa,2))/(2*pow(r,5))
                    *(1-(5*pow(curstates.k,2)/pow(r,2))); //constant value
        eci accgrav;
            accgrav.i=C*curstates.i;
            accgrav.j=C*curstates.j;
            accgrav.k=C*curstates.k;
        posstates dotstates;
            dotstates.i=curstates.vi+control.i; //ensure that these are ECI not ECEF
            dotstates.j=curstates.vj+control.j;
            dotstates.k=curstates.vk+control.k;
            dotstates.vi=-planet.sgp*curstates.i/pow(r,3)+accgrav.i+control.vi;
            dotstates.vj=-planet.sgp*curstates.j/pow(r,3)+accgrav.j+control.vj;
            dotstates.vk=-planet.sgp*curstates.k/pow(r,3)+accgrav.k+control.vk;
    };

}