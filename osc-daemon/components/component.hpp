#ifndef COMPONENT_H
#define COMPONENT_H

#include <string>
#include <iostream>

#include "../osctypes.hpp"

namespace osc {

  /** \class component
  Component class of a satellite object, used for the general definition
  of components on a satellite
  */
  class component{

    private:
    /// @param mass Mass of the object
    double mass;
    /// @param moi Moment of inertia of the object
    momentofinertia moi;
    /// @param pos Position of the object
    position pos;
    /// @param rot Quaternion rotation of the object
    quaternion rot;
    /// @param power Powermodel used by object
    powermodel power;

    public:
    
    /// Constructors
    component(double initMass, momentofinertia initMoi, position initPos, quaternion initRot, powermodel initPower):
      mass(initMass), moi(initMoi), pos(initPos), rot(initRot), power(initPower) { }

    // Accessers
    /// returns \p mass
    double getMass() {
      return mass;
    }
    /// returns \p moi
    momentofinertia getMoi() {
      return moi;
    }
    
    /// returns \p pos
    position getPos() {
      return pos;
    }

    /// returns \p power
    powermodel getPower() {
      return power;
    }

    /** sets the mass of an object
    @param[in] argMass the mass the object is set to
    */
    void setMass(double argMass) {
      mass = argMass;
    }
    
    /** sets the position of an object
    @param[in] argPos the position the object is set to
    */
    void setPos(position argPos) {
      pos = argPos;
    }
  };
}

#endif // COMPONENT_H