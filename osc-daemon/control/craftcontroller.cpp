#include "craftcontroller.hpp"

namespace osc {

    bool craftcontroller::initModel() {
        std::string pathString;
        std::cout << "Enter path to craft configuration: ";
        std::cin >> pathString;

        // craftconfig config = parseJson(pathString);
        
        // if (!config.populated()) return false;
        return true;
    }

    double craftcontroller::getMaxThrust() {
        return maxThrust;
    }

    double craftcontroller::getTransferISP() {
        return transferISP;
    }

    void craftcontroller::recomputeComponentDeps() {
        /*
        Recomputes craft parameters from a change in component configuiration
        */
        cg = position(0,0,0);
        moi = momentofinertia();
        mass = 0;
        wetMass = 0;

        for (auto i=components.begin(); i!=components.end(); i++) {
            cg = ((i->second.getPos().multiply(i->second.getMass())) + \
                    (cg.multiply((mass + wetMass)))).divide \
                    (i->second.getMass() + mass + wetMass);
            mass += i->second.getMass();
        }
        
        for (auto i=components.begin(); i!=components.end(); i++) {
            component* thisComponent = &i->second; // Assign pointer to class for readability
            position offset = thisComponent->getPos() - cg;
            moi.addMass(thisComponent->getMoi(),thisComponent->getMass(),offset);
        }

        for (auto i=fueltanks.begin(); i!=fueltanks.end(); i++) {
            cg = ((i->second.getFuelPos().multiply(i->second.getFuelMass())) + \
                    (cg.multiply((mass + wetMass)))).divide \
                    (i->second.getFuelMass() + mass + wetMass);
            wetMass += i->second.getFuelMass();
        }

        for (auto i=thrusters.begin(); i!=thrusters.end(); i++) {
            if (i->second.getAttitudeControl()) {
                attitudeActuators.insert({i->first, ftModel(i->second.getThrustAxis()*i->second.getMaxThrust(), components[i->first].getPos())}); // Add ftModel and string identifier to attitudeActuators
            }
        }

        for (auto i=rotators.begin(); i!=rotators.end(); i++) {
            vec3 torqueVec = i->second.getDirection()*i->second.getTorque();
            attitudeActuators.insert({i->first, ftModel(0, 0, 0, torqueVec[0], torqueVec[1], torqueVec[2])}); // Add ftModel and string identifier to attitudeActuators
            std::cout << "_Debug_ rotator torque in x: " << torqueVec[0] << std::endl;
        }
    }

    std::map<std::string, double> craftcontroller::forcesToCommands(ftModel setpoint) {
        std::map<std::string, double> thrusterCommands;

        bool matched = false;
        ftModel sp = setpoint;
        int count = 0;

        while (!matched) {
            int currAxis = sp.getDominantAxis();
            std::string actuatorID;
            for (auto i=attitudeActuators.begin(); i!= attitudeActuators.end(); i++) {
                if (i->second.getDominantAxis() == currAxis) {
                    double command = sp[currAxis]/i->second[currAxis]; // Eliminate command error in dominant axis

                    if (thrusterCommands.count(i->first)>0) {
                        thrusterCommands[i->first] = thrusterCommands[i->first] + command;
                    }

                    else {
                        thrusterCommands.insert({i->first, command});
                    }

                    sp = sp - (i->second * command);
                    count++;
                    break;
                }
            }
            if (sp.mag()<0.01*setpoint.mag() || count > 100) {
                matched = true;
            }
        }

        return thrusterCommands;
    }

    craftcontroller::craftcontroller() {
        if (!initModel()) throw std::runtime_error("Craft model failed to initialise");
    }

    void craftcontroller::beginControl() {
        // scheduler* schedule = new scheduler();
        std::cout << "Beginning control" << std::endl;
    }
};