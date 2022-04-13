#include <iostream>
#include <vector>
#include <map>

#include "../components/component.hpp"
#include "../components/fueltank.cpp"
#include "../osctypes.hpp"

namespace osc {

class craftcontroller {

    private:

    // Class vars (stack)
    double mass;
    double wetMass;
    position cg;
    momentofinertia moi;
    orbparam orbit;
    orbrot rotation;

    std::map<std::string, component> components;
    std::map<std::string, fueltank> fueltanks;

    // Member functions
    bool initModel() {
        // jsonModel model;
        // model.read("Path to json file");
        // if (model.hasData()) {
        //     components = model.getComponents();
        //     fueltanks = model.getFueltanks();
            return true;
        // }
        // else {
        //     return false;
        // }
    }

    void recomputeComponentDeps() {
        /*
        Recomputes craft parameters from a change in component configuiration
        */
        cg = position();
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
    }
    
    public:

    // Initialiser
    craftcontroller() {
        if (!initModel()) throw std::runtime_error("Craft model failed to initialise");
    }


    // Getters and setters
    // position getCG() {
    //     return cg;
    // }
};
}