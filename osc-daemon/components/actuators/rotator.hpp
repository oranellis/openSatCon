#ifndef ROTATOR_H
#define ROTATOR_H

#include <array>

namespace osc {

  /** \class rotator
  Rotator Class to store rotating attitude control systems, like reaction wheels or magnetorquers
  */
  class rotator {
    private:
      /** magnetorquer specific parameters, zero for a reaction wheel
      @param maxDipoleMoment maxmimum dipole moment
      */
      double maxDipoleMoment;
      
      /** rotation wheel specific parameters, zero for a magnetorquer
      @param torque torque exerted by reaction wheel
      @param storedMomentum momentum currently stored by reaction wheel
      @param maxRotationSpeed maximum rotation speed of reaction wheel
      */
      double torque;
      double storedMomentum;
      double maxRotationSpeed; 

      /** shared parameters
      @param moi Moment of inertia of the rotator
      @param direction Direction of vector of the axis of rotation
      */
      momentofinertia moi;
      vec3 direction;

    public:
      /// Constructor
      rotator(double initMaxDipoleMoment, double initTorque, double initStoredMomentum, double initMaxRotationSpeed, 
        momentofinertia initMoi, std::array<double, 3> initDireciton):
        maxDipoleMoment(initMaxDipoleMoment), 
        torque(initTorque), 
        storedMomentum(initStoredMomentum),
        maxRotationSpeed(initMaxRotationSpeed),
        moi(initMoi),
        direction(initDireciton) { }

      /// returns \p maxDipoleMoment
      double getMaxDipoleMoment() {
        return maxDipoleMoment; 
      }
      /// returns \p torque
      double getTorque () {
        return torque;
      }
      /// returns \p storedMomentum
      double getStoredMomentum() {
        return storedMomentum;
      }
      /// returns \p maxRotationSpeed
      double getMaxRotationSpeed() {
        return maxRotationSpeed;
      }
      /// returns \p moi
      momentofinertia getMomentOfInertia() {
        return moi;
      }
      /// returns \p direction
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