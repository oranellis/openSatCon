class Component{
  private:
    double Mass, Position[3], Power;
  public:
    double GetMass() {
      return Mass;
	}
    double * GetPos() {
      return *Position;
	}
    double GetPower() {
      return Power;
	}
    double MomentOfInertia(){
	} //TBD

    void SetMass(double mass) {
      Mass = mass;
	}
    void * SetPos(double position) {
      *Position = *position;
	}
};
