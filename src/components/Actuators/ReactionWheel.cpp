class ReactionWheel : Actuator{
  private:
    double Torque, StoredMomentum, RotationSpeed;
  public:
    double GetTorque(){
      return Torque;}
    double GetStoredMomentum(){
      return StoredMomentum;}
    double GetRotationSpeed(){
      return RotationSpeed;}
      //Torque  = RotationSpeed * Something ()
};