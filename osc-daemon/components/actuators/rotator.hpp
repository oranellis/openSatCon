#ifndef ROTATOR_H
#define ROTATOR_H

#include <array>

namespace osc {

  class rotator {
    private:
      //magnetorquer parameters
      double maxDipoleMoment;
      
      //rotation wheel parameters
      double torque;
      double storedMomentum;
      double maxRotationSpeed; 

      //shared parameters
      momentofinertia moi;
      vec3 direction;

    public:
      // Constructor
      rotator(double initMaxDipoleMoment, double initTorque, double initStoredMomentum, double initMaxRotationSpeed, 
        momentofinertia initMoi, std::array<double, 3> initDireciton):
        maxDipoleMoment(initMaxDipoleMoment), 
        torque(initTorque), 
        storedMomentum(initStoredMomentum),
        maxRotationSpeed(initMaxRotationSpeed),
        moi(initMoi),
        direction(initDireciton) { }


      double getMaxDipoleMoment() {
        return maxDipoleMoment;
      }
      double getTorque () {
        return torque;
      }
      double getStoredMomentum() {
        return storedMomentum;
      }
      double getMaxRotationSpeed() {
        return maxRotationSpeed;
      }
      momentofinertia getMomentOfInertia() {
        return moi;
      }
      vec3 getDirection() {
        return direction;
      }
  };

        //Torque  = DipoleMoment(Am^2) * MagneticFieldStrength(T)
        //        = DipoleMoment(Am^2) * (EarthMagneticMoment(Tm^3) / SatelliteRadius(m)^3 * lambda)

        //EarthMagneticMoment = 7.8e15Tm^3 (different for other planets)

        //lambda = sqrt(1 + (3 * sin(pi/2 - MagneticLatitude)^2))

        //MagneticLatitude = lookup table (International Geomagnetic Reference Field (IGRF))

}

#endif // ROTATOR_H