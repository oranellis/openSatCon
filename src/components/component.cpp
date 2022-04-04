#include "component.hpp"

#include "../osctypes.cpp"

class component{

  private:
  
    double mass;
    position datumPos;
    powermodel electro;


  public:

  //Initialisers
  component(double mass, double* position) {
    mass = mass;
    for (int i=0;i<3;i++) {
    position[i] = position[i];
    }
  }

  // Accessors
  double getMass() {
    return mass;
	}

  position getPos() {
    return datumPos;
	}

  void setMass(double argMass) {
    mass = argMass;
	}
  
  void setPos(double* position) {
    for (int i=0; i<3;i++) {
      position[i] = position[i];
    }
	}
};
