#ifndef COMPONENT_H
#define COMPONENT_H

class Component{
  private:
    double Mass, Position[3], Power;
  public:
    Component() {}
    double GetMass() {}
    double * GetPos() {}
    double GetPower() {}

    void SetMass(double mass) {}
    void * SetPos(double position) {}
};

#endif