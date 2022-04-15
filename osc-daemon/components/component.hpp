#ifndef COMPONENT_H
#define COMPONENT_H

#include <string>

#include "../osctypes.hpp"

namespace osc {

  class component{
    private:
      double mass;
      momentofinertia moi;
      position pos;
      powermodel power;
      std::string id; // Unique string identifier for the component

    public:
      component(double initMass, position initPos, powermodel initPower) {}
      double getMass() {}
      momentofinertia getMoi() {}
      position getPos() {}
      powermodel getPower() {}

      void setMass(double mass) {}
      void setPos(position position) {}
  };
}

#endif