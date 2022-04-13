#include <iostream>
#include <vector>
#include <math.h>

#include "osctypes.hpp"
#include "planet.cpp"

namespace osc{

    states dynamicmodel(states curstates, states control, inertiamatrix inverseinertia){//someone make matrices inversible pls
    states dotstates;
        dotstates.o1=inverseinertia.J11*control.o1+inverseinertia.J12*control.o2+inverseinertia.J13*control.o3;
        dotstates.o2=inverseinertia.J21*control.o1+inverseinertia.J22*control.o2+inverseinertia.J23*control.o3;
        dotstates.o3=inverseinertia.J31*control.o1+inverseinertia.J32*control.o2+inverseinertia.J33*control.o3;
        dotstates.q1=0.5*curstates.q1;
        dotstates.q2=0.5*curstates.q2;
        dotstates.q3=0.5*curstates.q3;
    };
}