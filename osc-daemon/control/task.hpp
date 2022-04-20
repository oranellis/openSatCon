#ifndef TASK_H
#define TASK_H

#include <chrono>
#include <iostream>
#include <array>

#include "../osctypes.hpp"
#include "craftcontroller.hpp"
#include "../orbitalmechanics/axistransforms.hpp"

namespace osc {

    /** \enum taskType
    enum describing the type of task
    */
    enum taskType {
        manoeuvre,
        groundTrack,
        example,
        empty
    };


    /** \class task
    \brief task definition
    Task class to define actions for the craft
    */
    class task {

        private:
        /// @param priority integer indicating priority of task, with higher integers corresponding to higher priority
        int priority;
        /// @param thisTaskType the type of task that the current instance represents
        taskType thisTaskType;
        /// @param startTime start time of the action
        std::chrono::time_point<std::chrono::system_clock> startTime;
        /// @param actionDuration duration of the action of the task
        std::chrono::microseconds actionDuration;
        /// @param KOE KOE parameters of the orbit
        orbParam KOE;
        /// @param posvelECI the eci coordinate system position and velocity
        eci posvelECI;
        /// @param pointVector the pointing vector in eci coordinates
        eci pointVector;
        /// @param timeOffset double of offset time
        double timeOffset;

        public:

        /** task() initialises an empty task*/
        task() {
            thisTaskType = taskType::empty;
        }

        /** task(transferISP, maxThrust, deltaV, impulseKOE, startMass) generates a new task with the input parameters
        @param[in] transferISP transfer ISP for a task
        @param[in] maxThrust the maximum thrust of the task
        @param[in] deltaV the deltaV of the task
        @param[in] impulseKOE the impulse KOE of the task
        @param[in] startMass the starting mass of the craft
        */
        task(double transferISP, double maxThrust, vnb deltaV, orbParam impulseKOE, double startMass) {//takes in KOE at impulse burn point, has mean anomaly
            pcs posvelPCSimpulse = KOEtoPCS(impulseKOE); // impulse koe at burn centre, converted perifocal coordinate system
            eci posvelECIimpulse = PCStoECI(impulseKOE, posvelPCSimpulse); // to eci, current pos and vel at impulse burn before burn

            double exhaustVel = transferISP*9.81;
            priority=10;

            double normDeltaV = deltaV.vVNB.mag();
            eci intermediatePointingVec;
            intermediatePointingVec = VNBtoECI(posvelECIimpulse, deltaV);

            actionDuration = std::chrono::microseconds((int)(startMass * exhaustVel / maxThrust
                             *(1 - exp(-normDeltaV / exhaustVel))*1000000));

            //trueAnom(burnTime-actionDuration/2)
            

            //rotation code (similar idea, find longest rotation time and perform task at that time before burnTime-actionDuration/2)
        };

        task(taskType type) {
            if (type == taskType::example) {
                thisTaskType = taskType::example;
            }
        }

        /** getPointingDirection() finds the pointing vector at an arbitrary time offset from a known pointing vector*/
        vec3 getPointingDirection() { //orbParam KOE, eci posvelECI, eci pointVector, double timeOffset

            vec3 pointingVec;

            switch (thisTaskType) {
                case taskType::manoeuvre: {
                    //this function will find the pointing vector at an arbitrary time offset from a known pointing vector
                    //a function to create a new pointing vector for multiple input types is explained below

                    orbParam offsetKOE = KOE; //have this equal KOE for everything except truAnom, other elements do not change while in motion

                    double meanAngularMotion = sqrt(planet.sgp/pow3(KOE.sma)); // this is mean anomaly (radians) travelled by the craft per second
                    // not accurate when ecc>0 therefore must convert to true anomaly

                    double meanAnomalyKnown = KOE.trueToMean(); // we find the mean anomaly at the known point
                    
                    double meanAnomalyOffset = meanAngularMotion * timeOffset; //then calculate and add the offset
                    double meanAnom = meanAnomalyKnown + meanAnomalyOffset; 

                    offsetKOE.meanToTrue(meanAnom); //this function then gives us the more accurate true anomaly

                    pcs posvelPCSoffset = KOEtoPCS(offsetKOE); //intermediate transform from KOE to ECI
                    eci posvelECIoffset = PCStoECI(offsetKOE, posvelPCSoffset); //gives ECI position and velocity at the offset point
                    // these position and velocity values are used for the ECItoVNB transform matrix,
                    // and to calculate the new pointing angle at the offset point;

                    eci newPointVector;
                    // newPointVector = posvelECI.rIJK - targetPosECI); //calculate new pointing vector, this is placeholder
                    // note that in some cases the target position in the ECI frame may have moved, and may need 
                    // recalculated, such as ground positions, fixed in ECEF, but moving in ECI.

                    //if targeting LLA position
                        //LLA->ECEF->ECI
                        //
                        //ecef posECEF = LLAtoECEF(posLLA);
                        //eci posECI = ECEFtoECI(posECEF, siderealTime) <--this needs to be the sidereal time at the offset time

                    //if targeting ground ECEF position
                        //ECEF->ECI
                        //
                        //eci posECI = ECEFtoECI(posECEF, siderealTime)

                    //if targeting point in space
                        //already done, its already an ECI position

                    //if doing an mnvr
                        //use relative velocity instead of relative position
                        //make sure that the time offset is worked out beforehand in that case, using this equation:
                        //actionDuration = std::chrono::microseconds((int)(startMass * exhaustVel / controller->getMaxThrust()*(1 - exp(-normDeltaV / exhaustVel))*1000000));

                    vnb newPointVectorVNB = ECItoVNB(posvelECIoffset, newPointVector); //translates and returns the ECI vector into a VNB vector
                    pointingVec = newPointVectorVNB.vVNB;
                    break;
                }

                case taskType::example: {
                    pointingVec = vec3(1,0,0);
                    break;
                }
            }

            return pointingVec;
        }

        /** getStartTime() returns \p startTime
        */
        auto getStartTime() {
            return startTime;
        }

        /** getTaskType() returns \p thisTaskType
        */
        taskType getTaskType() {
            return thisTaskType;
        }
    };

}

#endif // TASK_H