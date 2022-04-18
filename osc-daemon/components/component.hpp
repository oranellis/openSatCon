#ifndef COMPONENT_H
#define COMPONENT_H

#include <string>
#include <iostream>

#include "../osctypes.hpp"

namespace osc {

  class component{

    private:
    double mass;
    momentofinertia moi;
    position pos;
    quaternion rot;
    powermodel power;

    public:
    
    // Constructors
    component(double initMass, momentofinertia initMoi, position initPos, quaternion initRot, powermodel initPower):
      mass(initMass), moi(initMoi), pos(initPos), rot(initRot), power(initPower) { }

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