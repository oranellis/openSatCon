#include <stdlib.h>
#include <string>

class massobj {
    
    // Instance variables
    string identifier = "undefined";
    double mass = 0;
    double craftPos[3];

    // Inits
    massobj(string initIdentifier, double initMass, double initCraftPos[3]) {
        identifier = initCraftPos;
        mass = initMass;
        craftPos[0] = initCraftPos[0];
    }
}