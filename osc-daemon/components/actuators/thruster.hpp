#ifndef THRUSTER_H
#define THRUSTER_H

#include <array>

#include "../../osctypes.hpp"

namespace osc {

  class thruster {

    private:

    // Class variable initialisers
    double maxThrust; // (N)
    double minThrustFraction; // Value from 0 to 1 representing the minimum configurable thrust fraction
    double thrustFraction;
    double specificImpulse; // (s)
    bool attitudeControl;
    vec3 thrustAxis;
    position thrustCentre;

    public:
    //Constructor
    thruster(double initMaxThrust, double initMinFrac, double initThrustFrac, double initSpecificImpulse, 
    bool initAttitudeControl, vec3 initThrustAxis, position initThrustCentre):
      maxThrust(initMaxThrust), minThrustFraction(initMinFrac), thrustFraction(initThrustFrac),
      specificImpulse(initSpecificImpulse), attitudeControl(initAttitudeControl), thrustAxis(initThrustAxis),
      thrustCentre(initThrustCentre) { }

    double getMaxThrust() {
      return maxThrust;
    }

    double getSpecificImpulse() {
      return specificImpulse;
    }

    bool getAttitudeControl() {
      return attitudeControl;
    }

    vec3 getThrustAxis() {
      return thrustAxis;
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