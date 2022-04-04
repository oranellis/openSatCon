#include <iostream>
#include <vector>

#include "components/component.hpp"
#include "osctypes.cpp"

class craftcontroller {

    private:

    // Class vars (stack)
    std::vector<component> components;
    double dryMass;
    double wetMass;
    position cg;
    orbparam orbit;
    orbrot rotation;

    // Member functions
    bool initModel() {
        // read in the information from json file
        // parse info into relevant data structures
        return true;
    }

    void recomputeComponentDeps() {
        cg = position();
        dryMass = 0;
        wetMass = 0;

        for (int i=0; i<components.size(); i++) {
            
        }
    }
    
    public:

    // Initialiser
    craftcontroller() {
        if (!initModel()) throw std::runtime_error("Craft model failed to initialise");
    }


    // Getters and setters
    position getCG() {
        return cg;
    }
};