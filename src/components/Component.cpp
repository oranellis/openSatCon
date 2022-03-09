#include "Component.hpp"

class Component{

  private:
  
    double Mass;
    double Position[3];
    double Power;


  public:

  //Initialisers
  Component(double mass, double* position) {
    Mass = mass;
    for (int i=0;i<3;i++) {
      Position[i] = position[i];
    }
  }

  // Accessors
  double GetMass() {
    return Mass;
	}

  double* GetPos() {
    return Position;
	}

  void SetMass(double mass) {
    Mass = mass;
	}
  double GetPower() {
    return Power;
	}
  
  void SetPos(double* position) {
    for (int i=0; i<3;i++) {
      Position[i] = position[i];
    }
	}
};
