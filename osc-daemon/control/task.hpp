#ifndef TASK_H
#define TASK_H

#include <chrono>
#include <iostream>
#include <array>

#include "../osctypes.hpp"
#include "../control/craftcontroller.hpp"
#include "../orbitalmechanics/axistransforms.cpp"

namespace osc {
    class task {

        private:

        int priority;
        std::chrono::microseconds actionDuration;
        vec3 pointingVec; // this needs to be constantly changing throughout the burn so needs to be dynamic. Maybe make a seperate manoevure object that can be called for the current pointing vector that also includes the endpoint

        vec3 getPointingDirection() {
            vec3 pointingVec;
            
            impulseKOE.meanAnom = 2 * M_PI - (sqrt(planet.sgp / pow3(impulseKOE.sma)) * (actionDuration.count()/1000000) / 2);
            orbParam burnStartKOE = meanToTrue(impulseKOE); // converts to true anomaly
            pcs posvelPCSburnStart = KOEtoPCS(burnStartKOE);
            eci posvelECIburnStart = PCStoECI(impulseKOE, posvelPCSburnStart);
            //ECI parallel to pointingVecECI -> VNB transform of parallelECI <-this is where to point at start of burn
            pointingVec = ECItoVNB(posvelECIburnStart, intermediatePointingVec);
        }

        public:
        //initialiser
        
        task(craftcontroller *controller, vnb deltaV, orbParam impulseKOE, double startMass) {//takes in KOE at impulse burn point, has mean anomaly
            pcs posvelPCSimpulse = KOEtoPCS(impulseKOE); // impulse koe at burn centre, converted perifocal coordinate system
            eci posvelECIimpulse = PCStoECI(impulseKOE, posvelPCSimpulse); // to eci, current pos and vel at impulse burn before burn

            double exhaustVel = controller->getTransferISP()*9.81;
            priority=10;

            double normDeltaV = deltaV.vVNB.mag();
            eci intermediatePointingVec;
            intermediatePointingVec = VNBtoECI(posvelECIimpulse, deltaV);

            actionDuration = std::chrono::microseconds((int)(startMass * exhaustVel / controller->getMaxThrust()
                             *(1 - exp(-normDeltaV / exhaustVel))*1000000));

            //trueAnom(burnTime-actionDuration/2)
            

            //rotation code (similar idea, find longest rotation time and perform task at that time before burnTime-actionDuration/2)
        };
        task(vec3 rotAng, double trueAnom, double duration) {

        }

        vnb offsetVector(orbParam KOE, eci posvelECI, eci pointVector, double timeOffset){
            //this function will find the pointing vector at an arbitrary time offset from a known pointing vector

            orbParam offsetKOE; //have this equal KOE for everything except truAnom

            double meanAngularMotion = sqrt(planet.sgp/pow3(KOE.sma));

            double a = sqrt(1-pow2(KOE.ecc))*sin(KOE.truAnom);  
            double b = 1 + KOE.ecc*cos(KOE.truAnom);

            double meanAnomalyKnown = KOE.trueToMean(); // this is the best that I can get for this equation
            
            double meanAnomalyOffset = meanAngularMotion * timeOffset;

            double meanAnom = meanAnomalyKnown + meanAnomalyOffset; //this can be made into a double if required

            offsetKOE.meanToTrue(meanAnom);

            pcs posvelPCSoffset = KOEtoPCS(offsetKOE); //intermediate transform from KOE to ECI
            eci posvelECIoffset = PCStoECI(offsetKOE, posvelPCSoffset); //gives ECI position and velocity at the offset point
            // these position and velocity values are used for the ECItoVNB transform matrix,
            // and to calculate the new pointing angle at the offset point;

            eci newPointVector;
                newPointVector = posvelECI.rIJK.operator-(targetPosECI); //calculate new pointing vector
            // note that in some cases the target position in the ECI frame may have moved, and may need 
            // recalculated, such as ground positions, fixed in ECEF, but moving in ECI.

            vnb newPointVectorVNB = ECItoVNB(posvelECIoffset, newPointVector);
            return newPointVectorVNB;
        };
    };

}

#endif // TASK_H