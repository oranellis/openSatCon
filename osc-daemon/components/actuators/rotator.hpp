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
      position pos;
      std::array<double, 3> direction;

    public:
      // Constructor
      rotator(double initMaxDipoleMoment, double initTorque, double initStoredMomentum, double initMaxRotationSpeed,
       momentofinertia initMoi, position initPos, std::array<double, 3> initDireciton):
        maxDipoleMoment(initMaxDipoleMoment), torque(initTorque), storedMomentum(initStoredMomentum),
        maxRotationSpeed(initMaxRotationSpeed), moi(initMoi), pos(initPos), direction(initDireciton) { }


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
      position getPosition() {
        return pos;
      }
      std::array<double, 3> getDirection() {
        return direction;
      }
  };

  class magneticTorquer {
    private:
      double MagneticFluxDensity;
    public:
      double GetMagneticFluxDensity(){
        return MagneticFluxDensity;}
        //Torque  = DipoleMoment(Am^2) * MagneticFieldStrength(T)
        //        = DipoleMoment(Am^2) * (EarthMagneticMoment(Tm^3) / SatelliteRadius(m)^3 * lambda)

        //EarthMagneticMoment = 7.8e15Tm^3 (different for other planets)

        //lambda = sqrt(1 + (3 * sin(pi/2 - MagneticLatitude)^2))

        //MagneticLatitude = lookup table (International Geomagnetic Reference Field (IGRF))
  };

  class reactionWheel {
    private:
      double Torque, StoredMomentum, RotationSpeed;
      std::array<double,3> normal;
    public:
      double GetTorque(){
        return Torque;}
      double GetStoredMomentum(){
        return StoredMomentum;}
      double GetRotationSpeed(){
        return RotationSpeed;}
        //Torque  = RotationSpeed * Something ()
  };

  class rotator {

  };
}

#endif // ROTATOR_H