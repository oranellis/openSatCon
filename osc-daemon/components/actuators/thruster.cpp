class thruster {
  private:
    double Thrust, SpecificImpulse;
    double ThrustAxis[3];

  public:
    double GetThrust() {
      return Thrust;
      }

    double GetSpecificImpulse() {
      return SpecificImpulse;
      }
  //Torque = Thrust(N) * MomentArm(m)
};