#include "component.hpp"

class fueltank : public component{
    double fuelMass;
  public:
    double getFuelMass(){
      return fuelMass;}
    double getMass(){
      return getMass()+fuelMass;} //broken
    void burnFuel(double burnedFuel){
      fuelMass = fuelMass-burnedFuel; //assuming no CoG change for fuel use as of now
    }
}; //changing masses