#ifndef CRAFTCONTROLLER_H
#define CRAFTCONTROLLER_H

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "../osctypes.hpp"
#include "../datainjest/jsonparser.hpp"
#include "../components/component.hpp"
#include "../components/fueltank.hpp"
#include "../components/actuators/thruster.hpp"
#include "../components/actuators/rotator.hpp"
#include "scheduler.cpp"

namespace osc {

    class craftcontroller {

        private:

        // Class vars (stack)
        double mass;
        double wetMass;
        double maxThrust;
        position cg;
        momentofinertia moi;
        orbparam orbit;
        quaternion rotation;

        std::map<std::string, component> components;
        std::map<std::string, fueltank> fueltanks;
        std::map<std::string, thruster> thrusters;
        std::map<std::string, rotator> rotators;

        // Member functions
        
        public:

        // Initialiser
        craftcontroller();

        bool initModel();

        void recomputeComponentDeps();

        // double nominalManoevureThrust();

        // double minRotationRate();
        
        void beginControl();
    };
}

#endif // CRAFTCONTROLLER_H