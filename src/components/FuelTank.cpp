class FuelTank : Component{
    double FuelMass;
  public:
    double GetFuelMass(){
      return FuelMass;}
    double GetMass(){
      return Component::Mass+FuelMass;} //broken
    void BurnFuel(double BurnedFuel){
      FuelMass = FuelMass-BurnedFuel; //assuming no CoG change for fuel use as of now
    }
}; //changing masses