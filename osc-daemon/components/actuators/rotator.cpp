#include <array>

namespace osc {

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