#include "component.hpp"

class component{

  private:

    componentType type;
    double mass;
    momentofinertia moi;
    position pos;
    powermodel power;


  public:

  //Initialisers
  component(double initMass, position initPos, powermodel initPower):mass(initMass), pos(initPos), power(initPower) {
  }

  componentType getType () {
    return type;
  }

  // Accessors
  double getMass() {
    return mass;
	}

  momentofinertia getI() {
    return moi;
  }

  position getPos() {
    return pos;
	}

  void setMass(double argMass) {
    mass = argMass;
	}
  
  void setPos(position argPos) {
    for (int i=0; i<3;i++) {
      pos = argPos;
    }
	}
};
