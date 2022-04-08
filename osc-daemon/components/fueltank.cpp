#include "../osctypes.hpp"

namespace osc {

  class fueltank {

    std::string fuelType;
    double fuelMass;
    double fuelCapacity;
    position fuelMassPos;

    public:

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