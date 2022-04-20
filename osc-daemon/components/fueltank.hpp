#ifndef SENSOR_H
#define SENSOR_H

#include "../osctypes.hpp"

namespace osc {

  /** \class fueltank
  \brief Class for fuel tanks
  Fuel tank class
  */
  class fueltank {

    /// @param fuelType type of fuel used
    std::string fuelType;
    /// @param fuelMass mass of fuel on craft
    double fuelMass;
    /// @param fuelCapacity fuel capacity of fueltank
    double fuelCapacity;
    /// @param fuelMassPos position of center of mass of fuel
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

    /// returns \p fuelType
    std::string getFuelType() {
      return fuelType;
    }
    /// returns \p fuelMass
    double getFuelMass() {
      return fuelMass;
    }
    /// returns \p fuelCapacity
    double getFuelCapacity() {
      return fuelCapacity;
    }
    /// returns \p fuelMassPos
    position getFuelPos() {
      return fuelMassPos;
    }
    /** function to burn fuel mass
    @param[in] burnedFuel amount of fuel being burned
    @param[out] fuelMass outputs the new fuelmass value
    */
    void burnFuel(double burnedFuel){
      fuelMass = fuelMass-burnedFuel; //assuming no CoG change for fuel use as of now
    }
  }; //changing masses
}

#endif // SENSOR_H