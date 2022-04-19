#include <iostream>
#include <vector>
#include <math.h>

#include "../osctypes.hpp"
#include "planet.cpp"
#include "axistransforms.cpp"
#include "../control/task.hpp"

namespace osc {


std::vector<task> highlevelcommand(orbParam curKOE, orbParam aftKOE) {
    //this command will take in a set of desired orbital parameters and attempt to automatically create an efficicient 
    //burn to reach them
    std::vector<task> taskList;
    
    if (curKOE.ecc == 0 && aftKOE.ecc == 0) {
        bool HoBE = circOrbitChoice(curKOE, aftKOE);
        if (HoBE == true) {
            std::vector<task> taskListHohmann = hohmannTransfer(curKOE, aftKOE);//first burn at periapsis for efficiency
            for (int i=0; i<taskListHohmann.size(); i++) { taskList.push_back(taskListHohmann.at(i)); } // Adds all elements in list to taskList
        }
        else {
            std::vector<task> taskListBielliptic = biellipticTransfer(curKOE, aftKOE);
            for (int i=0; i<taskListBielliptic.size(); i++) { taskList.push_back(taskListBielliptic.at(i)); } // Adds all elements in list to taskList
        }//first burn at periapsis for efficiency
    }
    else {
        std::vector<task> taskListHohmann = hohmannTransfer(curKOE, aftKOE);
        for (int i=0; i<taskListHohmann.size(); i++) { taskList.push_back(taskListHohmann.at(i)); } // Adds all elements in list to taskList
    }

    if (curKOE.inc != aftKOE.inc) {
        std::vector<task> taskPlaneChange = planeChangeTransfer(curKOE, aftKOE);//perform this at the largest possible radius (apoapsis)
        for (int i=0; i<taskPlaneChange.size(); i++) { taskList.push_back(taskPlaneChange.at(i)); } // Adds all elements in list to taskList
    };
    //make big task list and include plane change modification
    return taskList;
};

std::vector<task> hohmannTransfer(orbParam curKOE, orbParam aftKOE) {
//a Hohmann transfer is a series of two burns in order to change both the periapse and apoapsis of an orbit.
//it is both deltaV and time efficient
    std::vector<task> taskListHohmann;
    vnb burn1dV, burn2dV;

    double rp1, ra1, rp2, ra2;
    rp1 = curKOE.sma * (1 - curKOE.ecc);
    ra1 = curKOE.sma * (1 + curKOE.ecc);
    rp2 = aftKOE.sma * (1 - aftKOE.ecc);
    ra2 = aftKOE.sma * (1 + aftKOE.ecc);

    double  h1 = angularMomentum(rp1, ra1);
    double  h2 = angularMomentum(rp2, ra2);
    double  h3 = angularMomentum(rp1, ra2);
    double  h3_= angularMomentum(ra1, rp2);

    double va1, vb2, va_1, vb_2, va3, vb3, va_3_, vb_3_;
    va1   = h1  / rp1;
    vb2   = h2  / ra2;
    va_1  = h1  / ra1;
    vb_2  = h2  / rp2;
    va3   = h3  / rp1;
    vb3   = h3  / ra2;
    va_3_ = h3_ / ra1;
    vb_3_ = h3_ / rp2;

    double dVa, dVb, dVa_, dVb_;
    dVa  = va3   - va1; 
    dVb  = vb2   - vb3;
    dVa_ = va_3_ - va_1;
    dVb_ = vb_2  - vb_3_;

    double  dV = abs(dVa)  + abs(dVb);
    double  dV_= abs(dVa_) + abs(dVb_);
    if (dV <= dV_) {
        burn1dV.vVNB.data[0] = dVa;
        burn2dV.vVNB.data[0] = dVb;

        task burnOne = task(burn1dV, 0); // taskListHohmann.push_back(task(burn1dV, 0))
        task burnTwo = task(burn2dV, M_PI);
    }    //do first burn (dVa)  at the current periapsis, circularise (dVb) at intermediate apoapsis
    else {
        burn1dV.vVNB.data[0] = dVa_;
        burn2dV.vVNB.data[0] = dVb_;

        task burnOne = task(burn1dV, M_PI);
        task burnTwo = task(burn2dV, 0);

    }          //do first burn (dVa_) at the current apoapsis, circularise (dVb_) at intermediate apoapsis

    //make taskList
    return taskListHohmann;
};

std::array<task> biellipticTransfer(orbParam curKOE, orbParam aftKOE) {
//A bielleptic transfer consists of three burns, with one large inital burn, followed by two correctional burns
//while it can be more deltaV efficient, it is much slower than a Hohmann transfer
    std::vector<task> taskListBielliptic;
    vnb burn1dV, burn2dV, burn3dV;
    double r1, rp2, ra2, rp3, ra3, r4;

    r1  = curKOE.sma;//feel free to optimise this
    rp2 = r1;
    ra2 = 2 * curKOE.sma;
    rp3 = aftKOE.sma;
    ra3 = ra2;
    r4  = rp3;

    double  va1 = sqrt(planet.sgp / r1);
    double  h2  = angularMomentum(rp2, ra2);
    double  h3  = angularMomentum(rp3, ra3);
    double  vc4 = sqrt(planet.sgp / r4);

    double va2, vb2, vb3, vc3;
    va2 = h2 / rp2;
    vb2 = h2 / ra2;
    vb3 = h3 / ra3;
    vc3 = h3 / rp3;


    burn1dV.vVNB.data[0] = va2 - va1; //prograde    at periapsis (true longitude=0)
    burn2dV.vVNB.data[0] = vb3 - vb2; //prograde    at apoapsis
    burn3dV.vVNB.data[0] = vc4 - vc3; //retrograde  at periapsis

    task burnOne   = task(burn1dV, 0); // taskListBielliptic.push_back(task(burn1dV, 0))
    task burnTwo   = task(burn2dV, M_PI);
    task burnThree = task(burn3dV, 0);
    //taskList
    return taskListBielliptic;
};

    // double massburned(double dV, double mo, double Isp){
    //     double propuse;

    //     double mf = mo * exp(-dV/Isp);
    //     propuse=mo-mf;
    //     return propuse;
    // };//calculate propellant mass used for a given delta V

double angularMomentum(double rp, double ra) {
    //just a useful function to have
    double h = sqrt(2 * planet.sgp) * sqrt((ra * rp) / (ra + rp));
    return h;
}

bool circOrbitChoice(orbParam curKOE, orbParam aftKOE) {
    //above a ratio of orbit radii of 15.58, the Bielliptic transfer is always more deltaV efficient
    //between 11.94 and 15.58, the Bielliptic will only be more efficient if its first burn is above
    //a certain radius, which increases as the ratio of radii appoaches 11.94, below which the 
    //Hohmann transfer is always more efficent
    double rc = aftKOE.sma;//final circle
    //double rb; apoapsis of biellipse
    double ra = curKOE.sma;//initial circle
     
    if (rc / ra < 11.94) {
        return true;
    }//hohmann transfer is better in this case

    else if(rc / ra > 15.58) {
        return false;
    }//bielliptic is better in this case
    
    else {
        return true;
    };  //bielliptic can be more dV efficient, but takes far longer
    //this will require a user choice, then take it as a function argument
    //double a=rc/ra;
    //double b=rb/ra;
    //double dVh=(1/sqrt(a))-((M_SQRT1_2*(1-a))/sqrt(a*(1+a)))-1;
    //double dVbe=sqrt((2*(a+b))/(a*b))-((1+sqrt(a))/(sqrt(a)))-sqrt(2/(b*(1+b)))*(1-b);
};

std::array<task> phasingTransfer(orbParam curKOE, double phasePeriod) {
    //this is a series of two burns that seek to change the timing of the orbit,
    //without affecting the positioning of it by any more than is necessary.
    std::vector<task> taskListPhasing;
    orbParam aftKOE;
    vnb burndV;
    double deltaV;
    aftKOE.sma = pow(((phasePeriod * sqrt(planet.sgp)) / (2 * M_PI)), (2/3));

    double ra, rb, rc;

    ra = curKOE.sma * (1 - curKOE.ecc);   //initial periapsis
    rb = curKOE.sma * (1 + curKOE.ecc);   //initial apoapsis 
    rc = 2 * aftKOE.sma - ra;             //phasing apoapsis

    double h1 = angularMomentum(ra, rb);
    double h2 = angularMomentum(ra, rc);
    double va1, va2, dV;

    va1                 = h1 / ra;
    va2                 = h2 / ra;
    burndV.vVNB.data[0] = va2 - va1;//first burn at periapsis
    //burns are prograde/retrograde if positive/negative
    task burnOne   = task(burndV, 0);
    task burnTwo   = task(-burndV, 2*M_PI);

    return taskListPhasing;
};

std::array<task> planeChangeTransfer(orbParam curKOE, orbParam aftKOE) {
    //this burn allows for the inclination of the orbit to be changed.
    //it is most deltaV efficient to perform this burn at the lowest
    //orbital veloctiy
    std::vector<task> taskPlaneChange;
    vnb deltaV;
    double deltaInc = aftKOE.inc - curKOE.inc;
    double r        = (curKOE.sma * (1 - pow2(curKOE.ecc))) / (1 + curKOE.ecc * cos(curKOE.truAnom));
    double curV     = sqrt(planet.sgp * ((2 / r) - (1 / curKOE.sma)));

    deltaV.vVNB.data[0] = curV * sin(deltaInc);
    deltaV.vVNB.data[1] = curV * (1 - cos(deltaInc));
    deltaV.vVNB.data[2] = 0;
        
    task taskPlaneChange = task(deltaV, M_PI); //perform at min velocity for delta V efficiency

    return taskPlaneChange;
};

task complexManeuver(double dVv, double dVn, double dVb, double burnTrueAnom) {
//this allows the user to choose their own burns at any vector they desire, which
//greatly increases the flexibility of the system over only simple burns
    std::vector<task> taskPlaneChange;
    vnb deltaV;

    double burnAngle = burnTrueAnom;

    deltaV.vVNB.data[0] = dVv;
    deltaV.vVNB.data[1] = dVn;
    deltaV.vVNB.data[2] = dVb;

    return task(deltaV, burnAngle); //return task with VNB values and req truanom
// };
};
}