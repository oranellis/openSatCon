#include <iostream>
#include <vector>
#include <map>

#include "components/component.hpp"
#include "components/fueltank.cpp"
#include "osctypes.hpp"

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
        // read in the information from json file
        // parse info into relevant data structures
        return true;
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