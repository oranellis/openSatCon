#ifndef CRAFTCONTROLLER_H
#define CRAFTCONTROLLER_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <thread>
#include <wiringPi.h>

#include "../osctypes.hpp"
#include "scheduler.hpp"
#include "task.hpp"
#include "../datainjest/jsonparser.hpp"
#include "../components/component.hpp"
#include "../components/fueltank.hpp"
#include "../components/actuators/thruster.hpp"
#include "../components/actuators/rotator.hpp"

namespace osc {

    /** \class craftcontroller
    \brief controller object
    Craft Controller class that creates a controller object 
    */
    class craftcontroller {

        private:

        /// @param CONTROL_LOOP_FREQ frequency of control loop
        const int CONTROL_LOOP_FREQ = 8000;
        /// @param KP KP of control loop
        const double KP = 0.01;

        const int PIN_X_POS = 16;
        const int PIN_X_NEG = 26;
        const int PIN_Y_POS = 5;
        const int PIN_Y_NEG = 6;
        const int PIN_Z_POS = 13;
        const int PIN_Z_NEG = 19;
        const int PIN_T = 1;

        // /** \fn controlLoopThread(*interupt, *curTask, *controller)
        // @param[in] interupt boolean value whether to interrupt current task
        // @param[in] curTask current task
        // @param[in] controller reference to craft controller
        // PID control loop 
        // */

        // Class vars (stack)
        /// @param mass Overall mass of craft
        double mass;
        /// @param wetMass fuel mass available
        double wetMass;
        /// @param maxThrust maximum thrust output
        double maxThrust;
        /// @param transferISP required ISP for the orbital transfer
        double transferISP;
        /// @param cg position of center of gravity of craft
        position cg;
        /// @param moi moment of inertia of craft
        momentofinertia moi;
        /// @param orbit Description of current orbit
        orbParam orbit;
        /// @param rotation Description of current rotation
        quaternion rotation;
        /// @param currTask The currently running task
        task currTask;
        /// @param taskInterupt An interupt flag for the running task, accessible from threads
        bool taskInterupt;

        /// @param components Map of all component objects to associated IDs
        std::map<std::string, component> components;
        /// @param fueltanks Map of all fueltank objects to associated IDs
        std::map<std::string, fueltank> fueltanks;
        /// @param thrusters Map of all thruster objects to associated IDs
        std::map<std::string, thruster> thrusters;
        /// @param rotators Map of all rotator objects to associated IDs
        std::map<std::string, rotator> rotators;
        /// @param attitudeActuators Map of all ftModels to associated IDs
        std::map<std::string, ftModel> attitudeActuators;
        /// @param thrusterCommands The commands to be sent ot the thusters as a double between 0 and 1, mapped to associated IDs
        std::map<std::string, double> thrusterCommands;

        // Member functions
        
        public:
        // Member functions
        bool initModel();

        double getMaxThrust();

        double getTransferISP();

        quaternion getRotation();

        void setRotation(quaternion argRotation) { rotation = argRotation; }

        void recomputeComponentDeps();
    
        void forcesToCommands(ftModel setpoint);
        
        void controlLoopThread();

        void sensorThread();

        void outputThread();

        craftcontroller();

        void beginExampleControl();
    };
}

#endif // CRAFTCONTROLLER_H