#include <iostream>
#include <vector>
#include <math.h>

#include "osctypes.hpp"
#include "planet.cpp"
#include "axistransforms.cpp"

namespace osc::orbmnvrs{


void highlevelcommand(orbparam curKOE, orbparam aftKOE){
if (curKOE.ecc==0&&aftKOE.ecc==0){
    bool HoBE = circOrbitChoice(curKOE, aftKOE);
    if (HoBE == true)
        {double deltaV = hohmannTransfer(curKOE, aftKOE);//first burn at periapsis for efficiency
}   else
        {double deltaV = biellipticTransfer(curKOE, aftKOE);}//first burn at periapsis for efficiency
}
else{double deltaV = hohmannTransfer(curKOE, aftKOE);}

if (curKOE.inc!=aftKOE.inc){
    vnb deltaV = planeChangeTransfer(curKOE, aftKOE);//perform this at the largest possible radius (apoapsis)
};
};

double hohmannTransfer(orbparam curKOE, orbparam aftKOE){
    double MNVRdV;
    double rp1, ra1, rp2, ra2;
    rp1=curKOE.sma*(1-curKOE.ecc);
    ra1=curKOE.sma*(1+curKOE.ecc);
    rp2=aftKOE.sma*(1-aftKOE.ecc);
    ra2=aftKOE.sma*(1+aftKOE.ecc);
    double  h1 =angularMomentum(rp1, ra1);
    double  h2 =angularMomentum(rp2, ra2);
    double  h3 =angularMomentum(rp1, ra2);
    double  h3_=angularMomentum(ra1, rp2);

    double va1,vb2,va_1,vb_2,va3,vb3,va_3_,vb_3_;
    va1     =h1/rp1;
    vb2     =h2/ra2;
    va_1    =h1/ra1;
    vb_2    =h2/rp2;
    va3     =h3/rp1;
    vb3     =h3/ra2;
    va_3_   =h3_/ra1;
    vb_3_   =h3_/rp2;

    double dVa, dVb, dVa_, dVb_;
    dVa =abs(va3-va1); //maybe if statement these for prograde/retrograde
    dVb =abs(vb2-vb3);
    dVa_=abs(va_3_-va_1);
    dVb_=abs(vb_2-vb_3_);

    double  dV =dVa+dVb;
    double  dV_=dVa_+dVb_;
    if (dV<=dV_){MNVRdV=dV;}    //do first burn (dVa)  at the current periapsis, circularise (dVb) at intermediate apoapsis
    else {MNVRdV=dV_;}          //do first burn (dVa_) at the current apoapsis, circularise (dVb_) at intermediate apoapsis
    return MNVRdV;
};

double biellipticTransfer(orbparam curKOE, orbparam aftKOE){
    double MNVRdV;
    double r1, rp2, ra2, rp3, ra3, r4;
    r1=curKOE.sma;//feel free to optimise this
    rp2=r1;
    ra2=2*curKOE.sma;
    rp3=aftKOE.sma;
    ra3=ra2;
    r4=rp3;
    double  va1=sqrt(planet.sgp/r1);
    double  h2=angularMomentum(rp2, ra2);
    double  h3=angularMomentum(rp3, ra3);
    double  vc4=sqrt(planet.sgp/r4);

    double va2, vb2, vb3, vc3;
    va2 =h2/rp2;
    vb2 =h2/ra2;
    vb3 =h3/ra3;
    vc3 =h3/rp3;


    double dVa, dVb, dVc;
    dVa =abs(va2-va1);//prograde    at periapsis (true longitude=0)
    dVb =abs(vb3-vb2);//prograde    at apoapsis
    dVc =abs(vc4-vc3);//retrograde  at periapsis

    MNVRdV=dVa+dVb+dVc;
    return MNVRdV;
};

    // double massburned(double dV, double mo, double Isp){
    //     double propuse;

    //     double mf = mo * exp(-dV/Isp);
    //     propuse=mo-mf;
    //     return propuse;
    // };//calculate propellant mass used for a given delta V

double angularMomentum(double rp, double ra){
    double h=sqrt(2*planet.sgp)*sqrt((ra*rp)/(ra+rp));
    return h;
}

bool circOrbitChoice(orbparam curKOE, orbparam aftKOE){
    double rc=aftKOE.sma;//final circle
    //double rb; apoapsis of biellipse
    double ra=curKOE.sma;//initial circle
     
    if (rc/ra<11.94){return true;}//hohmann transfer is better in this case
    else if(rc/ra>15.58){return false;}//bielliptic is better in this case
    else{return true;};//bielliptic can be more dV efficient, but takes far longer
    //this will require a user choice
    //double a=rc/ra;
    //double b=rb/ra;
    //double dVh=(1/sqrt(a))-((M_SQRT1_2*(1-a))/sqrt(a*(1+a)))-1;
    //double dVbe=sqrt((2*(a+b))/(a*b))-((1+sqrt(a))/(sqrt(a)))-sqrt(2/(b*(1+b)))*(1-b);
};

double phasingTransfer(orbparam curKOE, double phaseperiod){
    orbparam aftKOE;
    double deltaV;
    aftKOE.sma=pow(((phaseperiod*sqrt(planet.sgp))/(2*M_PI)),(2/3));
    double ra, rb, rc;
    ra=curKOE.sma*(1-curKOE.ecc);   //initial periapsis
    rb=curKOE.sma*(1+curKOE.ecc);   //initial apoapsis 
    rc=2*aftKOE.sma-ra;             //phasing apoapsis
    double h1=angularMomentum(ra, rb);
    double h2=angularMomentum(ra, rc);
    double va1, va2, dV1, dV2;
    va1=h1/ra;
    va2=h2/ra;
    dV1=va2-va1;//first burn at periapsis
    dV2=-dV1;//second burn at periapsis
    //burns are prograde/retrograde if positive/negative
    deltaV=abs(2*dV1);
    return deltaV;
};

vnb planeChangeTransfer(orbparam curKOE, orbparam aftKOE){
    vnb deltaV;
    double deltainc=aftKOE.inc-curKOE.inc;
    double r = (curKOE.sma*(1-pow2(curKOE.ecc)))/(1+curKOE.ecc*cos(curKOE.truanom));
    double curV=sqrt(planet.sgp*((2/r)-(1/curKOE.sma)));
    deltaV.v=curV*sin(deltainc);
    deltaV.n=curV*(1-cos(deltainc));
    deltaV.b=0;
    return deltaV; //perform at apoapsis for delta V efficiency
};

vnb complexManeuver(double dVv, double dVn, double dVb, double theta_burn){
    vnb deltaV;
    double burnAngle=theta_burn;
    deltaV.v=dVv;
    deltaV.n=dVn;
    deltaV.b=dVb;
    return deltaV; //return task with VNB values and req truanom
};
};