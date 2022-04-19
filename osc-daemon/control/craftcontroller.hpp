#ifndef CRAFTCONTROLLER_H
#define CRAFTCONTROLLER_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <thread>

#include "../osctypes.hpp"
#include "scheduler.hpp"
#include "../datainjest/jsonparser.hpp"
#include "../components/component.hpp"
#include "../components/fueltank.hpp"
#include "../components/actuators/thruster.hpp"
#include "../components/actuators/rotator.hpp"

namespace osc {

    struct craftcontroller {

        private:

        // Class vars (stack)
        double mass;
        double wetMass;
        double maxThrust;
        double transferISP;
        position cg;
        momentofinertia moi;
        orbParam orbit;
        quaternion rotation;

        std::map<std::string, component> components;
        std::map<std::string, fueltank> fueltanks;
        std::map<std::string, thruster> thrusters;
        std::map<std::string, rotator> rotators;
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

        void beginControl();
    };
}

#endif // CRAFTCONTROLLER_H