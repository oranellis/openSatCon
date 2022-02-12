#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <complex>
//#include <QtMath>
//#include <QGenericMatrix>


int rad2deg(double deg){
    double rad;
    rad=180/M_PI;
    return rad;
}

int main ()
{
  double MomentOfInertia[3][3];
  int ReactionWheelConfigMatrix[3][4]={{1,0,0,1},{0,1,0,1},{0,0,1,1}};
  double omega_dot, TorqueCommand, TorqueDisturbances, TorqueActuators, ReactionWheelAngularMomentum;
}

