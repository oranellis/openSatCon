#ifndef COMPONENT_H
#define COMPONENT_H

class component{
  private:
    double mass, position[3], power;
  public:
    component() {}
    double getMass() {}
    double * getPos() {}
    double getPower() {}

    void setMass(double mass) {}
    void * setPos(double position) {}
};

#endif