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
            i->second.getMaxThrust();
        }
    }

    craftcontroller::craftcontroller() {
        if (!initModel()) throw std::runtime_error("Craft model failed to initialise");
    }

    void craftcontroller::beginControl() {
        // scheduler* schedule = new scheduler();
        std::cout << "Beginning control" << std::endl;
    }
};