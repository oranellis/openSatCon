#include "../osctypes.hpp"

class fueltank {

  double fuelMass;
  position fuelMassPos;

  public:
  double getFuelMass(){
    return fuelMass;
  }

  position getFuelPos() {
    return fuelMassPos;
  }

  void burnFuel(double burnedFuel){
    fuelMass = fuelMass-burnedFuel; //assuming no CoG change for fuel use as of now
  }
}; //changing masses