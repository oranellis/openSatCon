#include <iostream>
#include <vector>
#include <math.h>

#include "../osctypes.hpp"
#include "planet.cpp"

namespace osc{

    /** \fn rotationDynamicModel(curStates, control, inverseInertia)
    @param[in] curStates current state of the model
    @param[in] control control state of model
    @param[in] inverseInertia inverse of the model inertia matrix
    @param[out] dotStates dotStates of the rotation 
    Outputs a set of rotations
    */
    rotStates rotationDynamicModel(rotStates curStates, rotStates control, std::array<std::array<double, 3>, 3> inverseInertia) { //someone make matrices inversible pls
    //this model will handle the rotational dynamics of the craft, and has six state vectors: 3 body rates and three quaternion vector parts
        //a simple RK4 integrated can be used to save computational power while retaining accuracy
    rotStates dotStates;
        
        dotStates.omega.data[0] = inverseInertia[0][0]*control.omega.data[0]+inverseInertia[0][1]*control.omega.data[1]+inverseInertia[0][2]*control.omega.data[2];
        dotStates.omega.data[1] = inverseInertia[1][0]*control.omega.data[0]+inverseInertia[1][1]*control.omega.data[1]+inverseInertia[1][2]*control.omega.data[2];
        dotStates.omega.data[2] = inverseInertia[2][0]*control.omega.data[0]+inverseInertia[2][1]*control.omega.data[1]+inverseInertia[2][2]*control.omega.data[2];
        dotStates.q             = curStates.q.operator/(2);
    };

    posStates positiondynamicmodel(posStates curStates, double thrust, double Ve) { 
        //this is the orbit positing dynamic model, and allows for the positon of the craft ot be effectively estimated, 
        //again using an RK4 integrator. there are seven states, 3 positions, 3 velocities, and a mass state, as the burning
        //of propellant will affect the dynamics of the spacecraft. Attempts were made to calculate some constants a single 
        //time, as this model may require a very low runtime
    double r = curStates.r.mag();
    double v = curStates.v.mag();
    double C =  (3 * planet.J2 * planet.sgp * pow2(planet.sMa)) / (2 * pow(r, 5)); //constant values
    double D = 5 * pow2(curStates.r.data[2] / r);                                  //to decrease runtime
    double accRelConst = (-planet.sgp / pow(r, 1.5));

        vec3 accRel = curStates.r.operator*(accRelConst);

         vec3 accGrav;
            accGrav.data[0] = C * (D - 1) * curStates.r.data[0];
            accGrav.data[1] = C * (D - 1) * curStates.r.data[1];
            accGrav.data[2] = C * (D - 3) * curStates.r.data[2];

        vec3 accRest = accRel.operator+(accGrav); //acceleration of system with no disturbances

        vec3 accThrust = curStates.v.operator*(thrust / curStates.m / v);

        posStates dotStates;
            dotStates.r = curStates.v;
            dotStates.v = accRest.operator+(accThrust);
            dotStates.m = -thrust / Ve * curStates.m;
    };

}