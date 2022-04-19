#ifndef CRAFTCONTROLLER_H
#define CRAFTCONTROLLER_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <thread>

#include "../osctypes.hpp"
#include "scheduler.hpp"
#include "attitudecontroller.hpp"
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
        
        void craftcontroller::controlLoopThread(bool *interupt, task *curTask);

        void beginControl();
    };
}

#endif // CRAFTCONTROLLER_H