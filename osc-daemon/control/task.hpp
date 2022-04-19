#ifndef TASK_H
#define TASK_H

#include <chrono>
#include <iostream>
#include <array>

namespace osc {
    class task {

        private:

        int priority;
        std::chrono::microseconds actionDuration;
        std::array<double, 3>* pointingVec; // this needs to be constantly changing throughout the burn so needs to be dynamic. Maybe make a seperate manoevure object that can be called for the current pointing vector that also includes the endpoint

        public:
        //initialiser
        
        task(vnb deltaV, orbParam KOE) {//takes in KOE at impulse burn point
            priority=10;
            double normDeltaV = deltaV.vVNB.mag();
            pointingVec = VNBtoECI(unitVec(deltaV));
            actionDuration = thrust/mass
            //trueAnom(burnTime-actionDuration/2)
            //ECI parallel to pointingVecECI -> VNB transform of parallelECI <-this is where to point at start of burn

            //rotation code (similar idea, find longest rotation time and perform task at that time before burnTime-actionDuration/2)
        };
        task(vec3 rotAng, double trueAnom, double duration) {

        }
    };

}

#endif // TASK_H