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