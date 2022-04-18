#ifndef THRUSTER_H
#define THRUSTER_H

#include <array>

#include "../../osctypes.hpp"

namespace osc {

  struct forceTorqueModel {
    /*
    Represents the vector of forces and moments of an actuator at maximum actuation
    */
    
    // Member vars
    std::array<double, 6> forceTorqueVec {0, 0, 0, 0, 0, 0};
    
    // Initialisers
    forceTorqueModel(double Fx, double Fy, double Fz,double Txx, double Tyy, double Tzz):forceTorqueVec({Fx, Fy, Fz, Txx, Tyy, Tzz}) {
      /* Generic initialiser for the vector, directly assigning each component */
    }

    forceTorqueModel(double maxThrust, std::array<double, 3> thrustAxis, position cg) {
      /* Initialiser for thrusters with a thrust magnitude and direction of action */
      
    }

    // Implicit type converters
    operator std::array<double, 6> () const { 
      /* Adds support for implicit conversion from forceTorqueModel to std::array<double, 6> */
      return forceTorqueVec; 
    }

  };

  class thruster {
    private:
    double maxThrust; // (N)
    double minThrustFraction; // Value from 0 to 1 representing the minimum configurable thrust fraction
    double thrustFraction;
    double specificImpulse; // (s)
    bool attitudeControl;
    std::array<double, 3> thrustDirection;
    position thrustCentre;

    public:
    double getMaxThrust() {
      return maxThrust;
    }

    double getSpecificImpulse() {
      return specificImpulse;
    }

    bool getAttitudeControl() {
      return attitudeControl;
    }

    std::array<double, 3> getThrustDirection() {
      return thrustDirection;
    }

    position getThrustCentre() {
      return thrustCentre;
    }

    bool setThrustFraction(double argThrustFraction) {
      if (argThrustFraction<=1 && argThrustFraction>=minThrustFraction) {
        thrustFraction = argThrustFraction;
        return true;
      }
      else {
        return false;
      }
    }
  };
}

#endif // THRUSTER_H