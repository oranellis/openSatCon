#ifndef THRUSTER_H
#define THRUSTER_H

#include <array>

#include "../../osctypes.hpp"

namespace osc {

  /** \class thruster 
  Thruster class for thrust producing components
  */
  class thruster {

    private:
    
    // Class variable initialisers
    /// @param maxThrust maximum thrust in Newtons (N)
    double maxThrust; 
    /// @param minThrustFraction minimum thrust fraction, representing minimum configurable thrust, between 0 and 1
    double minThrustFraction; // Value from 0 to 1 representing the minimum configurable thrust fraction
    /// @param thrustFraction normal thrust fraction, between 0 and 1
    double thrustFraction;
    /// @param specificImpulse specific impulse of thruster (s)
    double specificImpulse; // (s)
    /// @param attitudeControl boolean to define whether the thruster is used for attitude control
    bool attitudeControl;
    /// @param thrustAxis Vector defining the axis of the thrust
    vec3 thrustAxis;
    /// @param thrustCentre position of the center of thrust
    position thrustCentre;

    public:
    /// Constructor
    thruster(double initMaxThrust, double initMinFrac, double initThrustFrac, double initSpecificImpulse, 
    bool initAttitudeControl, vec3 initThrustAxis, position initThrustCentre):
      maxThrust(initMaxThrust), minThrustFraction(initMinFrac), thrustFraction(initThrustFrac),
      specificImpulse(initSpecificImpulse), attitudeControl(initAttitudeControl), thrustAxis(initThrustAxis),
      thrustCentre(initThrustCentre) { }

    /// returns \p maxThrust
    double getMaxThrust() {
      return maxThrust;
    }
    /// returns \p specificImpulse
    double getSpecificImpulse() {
      return specificImpulse;
    }
    /// returns \p attitudeControl
    bool getAttitudeControl() {
      return attitudeControl;
    }
    /// returns \p thrustAxis
    vec3 getThrustAxis() {
      return thrustAxis;
    }
    /// returns \p thrustCentre
    position getThrustCentre() {
      return thrustCentre;
    }
    /// sets the thrust fraction to a defined value, returns true if the setting is possible
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