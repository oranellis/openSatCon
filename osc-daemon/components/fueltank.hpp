#ifndef SENSOR_H
#define SENSOR_H

#include "../osctypes.hpp"

namespace osc {

  class fueltank {

    std::string fuelType;
    double fuelMass;
    double fuelCapacity;
    position fuelMassPos;

    public:

    //constructors
    fueltank(std::string initfuelType, double initFuelMass, double initFuelCapacity, position initPosition):
      fuelType(initfuelType), fuelMass(initFuelMass), fuelCapacity(initFuelCapacity), fuelMassPos(initPosition) { }

    fueltank(std::string initFuelType, double initFuelMass, position initFuelMassPos) {
      fuelType = initFuelType;
      fuelMass = initFuelMass;
      fuelMassPos = initFuelMassPos;
    };

    std::string getFuelType() {
      return fuelType;
    }

    double getFuelMass() {
      return fuelMass;
    }

    double getFuelCapacity() {
      return fuelCapacity;
    }

    position getFuelPos() {
      return fuelMassPos;
    }

    void burnFuel(double burnedFuel){
      fuelMass = fuelMass-burnedFuel; //assuming no CoG change for fuel use as of now
    }
  }; //changing masses
}

#endif // SENSOR_H