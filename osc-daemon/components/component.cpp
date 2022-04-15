#include "component.hpp"

namespace osc {

  //Initialisers
  component::component(double initMass, position initPos, powermodel initPower):mass(initMass), pos(initPos), power(initPower) { }

  // Accessors
  double component::getMass() {
    return mass;
  }

  momentofinertia component::getMoi() {
    return moi;
  }

  position component::getPos() {
    return pos;
  }

  void component::setMass(double argMass) {
    mass = argMass;
  }
  
  void component::setPos(position argPos) {
    for (int i=0; i<3;i++) {
      pos = argPos;
    }
  }
} // namespace osc
