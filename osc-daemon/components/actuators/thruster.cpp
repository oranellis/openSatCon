#include <array>

#include "../../osctypes.hpp"

namespace osc {
  class thruster {
    private:
      double maxThrust; // (N)
      double minThrustFraction; // Value from 0 to 1 representing the minimum configurable thrust fraction
      double thrustFraction;
      double specificImpulse; // (s)
      std::array<double, 3> thrustAxis;

    public:
      double getMaxThrust() {
        return maxThrust;
        }

      double GetSpecificImpulse() {
        return specificImpulse;
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
    //Torque = Thrust(N) * MomentArm(m)
  };
}