#ifndef COMPONENT_H
#define COMPONENT_H

#include <string>
#include <iostream>

#include "../osctypes.hpp"

namespace osc {

  class component{

    private:
    std::string id; // Unique string identifier for the component
    double mass;
    momentofinertia moi;
    position pos;
    quaternion rot;
    powermodel power;

    public:
    
    // Constructors
    component(double initMass, position initPos, powermodel initPower):mass(initMass), pos(initPos), power(initPower) { }

    // Accessers
    double getMass() {
      return mass;
    }

    momentofinertia getMoi() {
      return moi;
    }
    
    position getPos() {
      return pos;
    }

    powermodel getPower() {
      return power;
    }

    void setMass(double argMass) {
      mass = argMass;
    }

    void setPos(position argPos) {
      pos = argPos;
    }
  };
}

#endif // COMPONENT_H