#include <iostream>
#include <vector>
#include <math.h>
#include <chrono>

#include "osctypes.hpp"
#include "planet.cpp"

namespace osc{

double meanToTrue(orbParam KOE){
    // converts mean anomaly to true anomaly
        if (KOE.ecc < 0.2){ //this is a handy calculation to save time for small eccentricities
            KOE.truAnom = KOE.meanAnom+2*KOE.ecc*sin(KOE.meanAnom)+1.25*pow2(KOE.ecc)*sin(2*KOE.meanAnom)-pow3(KOE.ecc)*(0.25*sin(KOE.meanAnom)-(13/12)*sin(3*KOE.meanAnom));
        } else if(KOE.ecc < 1.0){// newton raphson's method must be used for higher eccentricities, e>1 is a parabolic orbit
            KOE.eccAnom=KOE.meanAnom+((KOE.ecc*sin(KOE.meanAnom)/(cos(KOE.ecc)-(M_PI_2-KOE.ecc)*sin(KOE.ecc)+KOE.meanAnom*sin(KOE.ecc))));
            double dE=KOE.eccAnom;
            while(abs(dE)>10e-10){
                dE=(KOE.eccAnom-KOE.ecc*sin(KOE.eccAnom)-KOE.meanAnom)/(1-KOE.ecc*cos(KOE.eccAnom));
                KOE.eccAnom=KOE.eccAnom-dE;
            };
            KOE.truAnom=atan2(sqrt(1-pow2(KOE.ecc)*sin(KOE.ecc)),cos(KOE.eccAnom)-KOE.ecc);
        };
};


double greenwichSiderealAngle(){
    //calculates the angle required for the transformation between ECI and ECEF
    double gsraret;
    time_t timer;
    //for Oran
    struct tm J2000 = {0};
    double seconds;

    J2000.tm_hour = 11;  J2000.tm_min = 58; J2000.tm_sec = 55;
    J2000.tm_year = 100; J2000.tm_mon = 0; J2000.tm_mday = 1;

    time(&timer);

    seconds = difftime(timer,mktime(&J2000));
    double julianDays=seconds/60/60/365.25;
    gsraret=2*M_PI*(0.77905722732640+1.00273781191135448*(julianDays-2451545));
    return gsraret;
    //this could be made down to sub-second accuracy if possible
};


ecef LLAtoECEF(lla posLLA) {
    ecef posECEF;
    // converts the angular position of the satellite to the ECEF co-ordinate system
    // this gives the correct altitude of a non-spherical earth, note that altitude here is above the ground
    double normalDistance=planet.sMa/sqrt(1-(planet.ecc*sin(posLLA.lat)));
    posECEF.rXYZ.data[0]=(normalDistance+posLLA.alt)*cos(posLLA.lon)*cos(posLLA.lat);
    posECEF.rXYZ.data[1]=(normalDistance+posLLA.alt)*sin(posLLA.lon)*cos(posLLA.lat);
    posECEF.rXYZ.data[2]=(normalDistance*(1-pow2(planet.ecc))+posLLA.alt)*sin(posLLA.lat);
    return posECEF;
};

lla ECEFtoLLA(ecef posECEF){
    //returns ground position of sub satellite point and satellite altitude from ECEF co-ords
    lla llaret;
    double secondeccentricity=planet.ecc/sqrt(1-pow2(planet.ecc));
    double p = sqrt(pow2(posECEF.rXYZ.data[2])+pow2(posECEF.rXYZ.data[1]));
    double theta = atan2(posECEF.rXYZ.data[2]*planet.sMa, p*planet.sma);
        llaret.lon = atan2(posECEF.rXYZ.data[1], posECEF.rXYZ.data[0]);
        llaret.lat = atan2(posECEF.rXYZ.data[2]+pow2(secondeccentricity)*planet.sma*pow3(sin(theta)),p-pow2(planet.ecc)*planet.sMa*pow3(sin(theta)));
    double normaldistance = planet.sMa/sqrt(1-pow2(planet.ecc)*pow2(sin(llaret.lat)));
        llaret.alt = (p/cos(llaret.lat))-normaldistance;
    return llaret;
};

ned ECEFtoNED(ecef posECEF, lla refLLA){
    // returns vector from satellite to ground station in North East Down reference frame
    ned posNED;
    ecef refECEF = LLAtoECEF(refLLA);
    ecef relECEF;
    relECEF.rXYZ=posECEF.rXYZ.operator-(refECEF.rXYZ);

    posNED.rNED.data[0] =   -sin(refLLA.lat)*cos(refLLA.lon)*relECEF.rXYZ.data[0]+
                            sin(refLLA.lat)*sin(refLLA.lon)*relECEF.rXYZ.data[1]+
                            cos(refLLA.lat)*relECEF.rXYZ.data[2];

    posNED.rNED.data[1] =   -sin(refLLA.lon)*relECEF.rXYZ.data[0]+
                            cos(refLLA.lon)*relECEF.rXYZ.data[1]+
                            0;

    posNED.rNED.data[2] =   -cos(refLLA.lat)*cos(refLLA.lon)*relECEF.rXYZ.data[0]+
                            -cos(refLLA.lat)*sin(refLLA.lon)*relECEF.rXYZ.data[1]+
                            -sin(refLLA.lat)*relECEF.rXYZ.data[2];

    return posNED;
};

ecef NEDtoECEF(ned posNED, lla refLLA){
    // returns ECEF position from NED vector
    ecef posECEF;

    posECEF.rXYZ.data[0] = -cos(refLLA.lon)*sin(refLLA.lat)*posNED.rNED.data[0]
                           -sin(refLLA.lon)*posNED.rNED.data[1]
                           -cos(refLLA.lon)*cos(refLLA.lat)*posNED.rNED.data[2];

    posECEF.rXYZ.data[1] = -sin(refLLA.lon)*sin(refLLA.lat)*posNED.rNED.data[0]
                           +cos(refLLA.lon)*posNED.rNED.data[1]
                           -sin(refLLA.lon)*cos(refLLA.lat)*posNED.rNED.data[2];

    posECEF.rXYZ.data[2] = cos(refLLA.lat)*posNED.rNED.data[0]
                           +0
                           -sin(refLLA.lat)*posNED.rNED.data[2];
    return posECEF;   
};

ecef EARtoECEF(ear posEAR, lla refLLA){
    // returns ECEF position from ground station Elevation Azimuth Range tracking
    // this method goes through SEZ co-ordinates, but seems to be the best method
    ecef posECEF;
    ecef refECEF = LLAtoECEF(refLLA);
    posECEF.rXYZ.data[0] = sin(refLLA.lat)*cos(refLLA.lon)*(-posEAR.r*cos(posEAR.e)*cos(posEAR.a))
                           -sin(refLLA.lon)*(posEAR.r*cos(posEAR.e)*sin(posEAR.a))
                           +cos(refLLA.lat)*cos(refLLA.lon)*(posEAR.r*sin(posEAR.e)) + refECEF.rXYZ.data[0];
    
    posECEF.rXYZ.data[1] = sin(refLLA.lat)*sin(refLLA.lon)*(-posEAR.r*cos(posEAR.e)*cos(posEAR.a))
                           +cos(refLLA.lon)*(posEAR.r*cos(posEAR.e)*sin(posEAR.a))
                           +cos(refLLA.lat)*sin(refLLA.lon)*(posEAR.r*sin(posEAR.e)) + refECEF.rXYZ.data[1];

    posECEF.rXYZ.data[2] = -cos(refLLA.lat)*(-posEAR.r*cos(posEAR.e)*cos(posEAR.a))
                           +0
                           +sin(refLLA.lat)*(posEAR.r*sin(posEAR.e)) + refECEF.rXYZ.data[2];

    return posECEF;
};


ear ECEFtoEAR(ecef posECEF, lla refLLA){
    // returns the EAR values seen by a ground station of a satellite at the given ECEF position
    ear posEAR;
    ecef refECEF= LLAtoECEF(refLLA);
    ecef relECEF;
    relECEF.rXYZ=posECEF.rXYZ.operator-(refECEF.rXYZ);

    posEAR.e=   acos((refECEF.rXYZ.data[0]*relECEF.rXYZ.data[0]+refECEF.rXYZ.data[1]*relECEF.rXYZ.data[1]+refECEF.rXYZ.data[2]*relECEF.rXYZ.data[2])
                /sqrt((pow2(refECEF.rXYZ.data[0])+pow2(refECEF.rXYZ.data[1])+pow2(refECEF.rXYZ.data[2]))*(pow2(relECEF.rXYZ.data[0])+pow2(relECEF.rXYZ.data[1])+pow2(relECEF.rXYZ.data[2]))))-M_PI_2;

    posEAR.a=   atan((-refECEF.rXYZ.data[1]*relECEF.rXYZ.data[0]+refECEF.rXYZ.data[0]*relECEF.rXYZ.data[1]/sqrt((pow2(refECEF.rXYZ.data[0])+pow2(refECEF.rXYZ.data[1])+pow2(refECEF.rXYZ.data[2]))*(pow2(relECEF.rXYZ.data[0])+pow2(relECEF.rXYZ.data[1])+pow2(relECEF.rXYZ.data[2]))))
                /(-refECEF.rXYZ.data[2]*refECEF.rXYZ.data[0]*relECEF.rXYZ.data[0]-refECEF.rXYZ.data[2]*refECEF.rXYZ.data[1]*relECEF.rXYZ.data[1]+(pow2(refECEF.rXYZ.data[0])+pow2(refECEF.rXYZ.data[1]))*relECEF.rXYZ.data[2])/sqrt((pow2(refECEF.rXYZ.data[0])+pow2(refECEF.rXYZ.data[1]))+(pow2(refECEF.rXYZ.data[0])+pow2(refECEF.rXYZ.data[1])+pow2(refECEF.rXYZ.data[2]))*(pow2(relECEF.rXYZ.data[0])+pow2(relECEF.rXYZ.data[1])+pow2(relECEF.rXYZ.data[2]))));

    posEAR.r=   sqrt(pow2(relECEF.rXYZ.data[0])+pow2(relECEF.rXYZ.data[1])+pow2(relECEF.rXYZ.data[2]));
    return posEAR;
};

// thcs EARtoSEZ(ear argpos, ear argvel, eci refpos) {// need to double check this maths
//     //also probably better suited elsewhere as this is a tracking calculation
//     thcs sezposret;
//     thcs sezvelret;
//     sezposret.s=refpos.rIJK.data[0]- argpos.r*cos(argpos.e)*cos(argpos.a);
//     sezposret.e=refpos.rIJK.data[1]- argpos.r*cos(argpos.e)*sin(argpos.a);
//     sezposret.rXYZ.data[2]=refpos.rIJK.data[2]- argpos.r*sin(argpos.e);

//     sezvelret.s=-argvel.r*cos(argpos.e)*cos(argpos.a)+argpos.r*sin(argpos.e)*cos(argpos.a)*argvel.e+argpos.r*cos(argpos.e)*sin(argpos.a)*argvel.a;
//     sezvelret.e=argvel.r*cos(argpos.e)*sin(argpos.a)-argpos.r*sin(argpos.e)*sin(argpos.a)*argvel.e+argpos.r*cos(argpos.e)*cos(argpos.a)*argvel.a;
//     sezvelret.rXYZ.data[2]=argvel.r*sin(argpos.e)*argpos.r*cos(argpos.e)*argvel.e;
//     return sezposret;
//     return sezvelret;
// };

ecef ENUtoECEF(enu satpos, lla reflla){
    // returns the ECEF position of a satellite tracked by a ground station in East North Up
    ecef ecefret;
    ecef refpos = LLAtoECEF(reflla);

    ecefret.rXYZ.data[0]=  -sin(reflla.lon)*satpos.e
                -sin(reflla.lat)*cos(reflla.lon)*satpos.n
                +cos(reflla.lat)*cos(reflla.lon)*satpos.u + refpos.rXYZ.data[0];

    ecefret.rXYZ.data[1]=  cos(reflla.lon)*satpos.e
                -sin(reflla.lat)*sin(reflla.lon)*satpos.n
                +cos(reflla.lat)*sin(reflla.lon)*satpos.u + refpos.rXYZ.data[1];

    ecefret.rXYZ.data[2]=  0
                +cos(reflla.lat)*satpos.n
                +sin(reflla.lat)*satpos.u + refpos.rXYZ.data[2];

    return ecefret;
};

enu ECEFtoENU(ecef satpos, lla reflla){
    // returns the ENU values seen by a ground station of a satellite at the given ECEF position
    enu enuret;
    ecef refpos = LLAtoECEF(reflla);

    enuret.e=   -sin(reflla.lon)*(satpos.rXYZ.data[0]-refpos.rXYZ.data[0])
                +cos(reflla.lon)*(satpos.rXYZ.data[1]-refpos.rXYZ.data[1])
                +0;

    enuret.n=   -sin(reflla.lat)*cos(reflla.lon)*(satpos.rXYZ.data[0]-refpos.rXYZ.data[0])
                -sin(reflla.lat)*sin(reflla.lon)*(satpos.rXYZ.data[1]-refpos.rXYZ.data[1])
                +cos(reflla.lat)*(satpos.rXYZ.data[2]-refpos.rXYZ.data[2]);

    enuret.u=   cos(reflla.lat)*cos(reflla.lon)+(satpos.rXYZ.data[0]-refpos.rXYZ.data[0])
                +cos(reflla.lat)*sin(reflla.lon)*(satpos.rXYZ.data[1]-refpos.rXYZ.data[1])
                +sin(reflla.lat)*(satpos.rXYZ.data[2]-refpos.rXYZ.data[2]);

    return enuret;
};

ecef SEZtoECEF(thcs satpos, lla reflla){
    // returns the ECEF position of a satellite tracked by a ground station in the Topocentric Horizon Co-ordinate System
    ecef ecefret;
    ecef refpos = LLAtoECEF(reflla);

    ecefret.rXYZ.data[0]=  cos(reflla.lon)*satpos.s
                +sin(reflla.lat)*cos(reflla.lon)*satpos.e
                + refpos.rXYZ.data[2];

    ecefret.rXYZ.data[1]=  -sin(reflla.lat)*sin(reflla.lon)*satpos.s
                +sin(reflla.lat)*cos(reflla.lon)*satpos.e
                +cos(reflla.lat)*satpos.rXYZ.data[2] + refpos.rXYZ.data[1];

    ecefret.rXYZ.data[2]=  cos(reflla.lat)*sin(reflla.lon)*satpos.s
                -cos(reflla.lat)*cos(reflla.lon)*satpos.e
                +sin(reflla.lat)*satpos.rXYZ.data[2] + refpos.rXYZ.data[2];

    return ecefret;
};

thcs ECEFtoSEZ(ecef satpos, lla reflla){
    // returns the SEZ values seen by a ground station of a satellite at the given ECEF position
    thcs sezret;
    ecef refpos = LLAtoECEF(reflla);

    sezret.s=   sin(reflla.lat)*cos(reflla.lon)*(satpos.rXYZ.data[0]-refpos.rXYZ.data[0])
                +sin(reflla.lat)*sin(reflla.lon)*(satpos.rXYZ.data[1]-refpos.rXYZ.data[1])
                -cos(reflla.lat)*(satpos.rXYZ.data[2]-refpos.rXYZ.data[2]);

    sezret.e=   -sin(reflla.lon)*(satpos.rXYZ.data[0]-refpos.rXYZ.data[0])
                +cos(reflla.lon)*(satpos.rXYZ.data[1]-refpos.rXYZ.data[1])
                +0;

    sezret.rXYZ.data[2]=   cos(reflla.lat)*cos(reflla.lon)+(satpos.rXYZ.data[0]-refpos.rXYZ.data[0])
                +cos(reflla.lat)*sin(reflla.lon)*(satpos.rXYZ.data[1]-refpos.rXYZ.data[1])
                +sin(reflla.lat)*(satpos.rXYZ.data[2]-refpos.rXYZ.data[2]);

    return sezret;
};

eci PCStoECI(orbParam arg, pcs pcsarg){
    eci eciret;

    eciret.rIJK.data[0]=   (cos(arg.asc)*cos(arg.aop)-sin(arg.asc)*sin(arg.aop)*cos(arg.rIJK.data[0]nc))*pcsarg.p +
                (-cos(arg.asc)*sin(arg.aop)-cos(arg.asc)*sin(arg.aop)*cos(arg.rIJK.data[0]nc))*pcsarg.q +
                (sin(arg.asc)*sin(arg.rIJK.data[0]nc))*pcsarg.w;

    eciret.rIJK.data[1]=   (sin(arg.asc)*cos(arg.aop)-cos(arg.asc)*sin(arg.aop)*cos(arg.rIJK.data[0]nc))*pcsarg.p +
                (-sin(arg.asc)*sin(arg.aop)-cos(arg.asc)*cos(arg.aop)*cos(arg.rIJK.data[0]nc))*pcsarg.q +
                (-cos(arg.asc)*sin(arg.rIJK.data[0]nc))*pcsarg.w;

    eciret.rIJK.data[2]=   sin(arg.aop)*sin(arg.rIJK.data[0]nc)*pcsarg.p +
                cos(arg.aop)*sin(arg.rIJK.data[0]nc)*pcsarg.q +
                cos(arg.rIJK.data[0]nc)*pcsarg.w;

    eciret.vi=   (cos(arg.asc)*cos(arg.aop)-sin(arg.asc)*sin(arg.aop)*cos(arg.rIJK.data[0]nc))*pcsarg.vp +
                (-cos(arg.asc)*sin(arg.aop)-cos(arg.asc)*sin(arg.aop)*cos(arg.rIJK.data[0]nc))*pcsarg.vq +
                (sin(arg.asc)*sin(arg.rIJK.data[0]nc))*pcsarg.vw;

    eciret.vj=   (sin(arg.asc)*cos(arg.aop)-cos(arg.asc)*sin(arg.aop)*cos(arg.rIJK.data[0]nc))*pcsarg.vp +
                (-sin(arg.asc)*sin(arg.aop)-cos(arg.asc)*cos(arg.aop)*cos(arg.rIJK.data[0]nc))*pcsarg.vq +
                (-cos(arg.asc)*sin(arg.rIJK.data[0]nc))*pcsarg.vw;

    eciret.vk=   sin(arg.aop)*sin(arg.rIJK.data[0]nc)*pcsarg.vp +
                cos(arg.aop)*sin(arg.rIJK.data[0]nc)*pcsarg.vq +
                cos(arg.rIJK.data[0]nc)*pcsarg.vw;

    return eciret;
};

pcs KOEtoPCS(orbParam arg){
    // gives the ECI position and velocity of a satellite given its Kepler Orbital Elements

    pcs pcsret;
    double orbitradius;
    orbitradius = (arg.sma*(1-pow2(arg.ecc)))/(1+arg.ecc*cos(arg.truAnom));
    pcsret.p = orbitradius * cos(arg.truAnom); 
    pcsret.q = orbitradius * sin(arg.truAnom);
    pcsret.w = 0;

    double p = arg.sma*(1-pow2(arg.ecc));
    double thetadot = sqrt(planet.sgp*p)/pow2(orbitradius);
    double rdot = sqrt(planet.sgp/p)*arg.ecc*sin(arg.truAnom);

    pcsret.vp = rdot*cos(arg.truAnom)-orbitradius*sin(arg.truAnom)*thetadot;
    pcsret.vq = rdot*sin(arg.truAnom)-orbitradius*cos(arg.truAnom)*thetadot;
    pcsret.vw = 0;
    return pcsret;
};

orbParam ECItoKOE(eci posvel){
    orbParam koeret;
    eci argangmnt;
    eci nodalvect;
    pcs eccvect;
    double orbitalradius = sqrt(pow2(posvel.rIJK.data[0])+pow2(posvel.rIJK.data[1])+pow2(posvel.rIJK.data[2]));
    double orbitalvelocity = sqrt(pow2(posvel.vi)+pow2(posvel.vj)+pow2(posvel.vk));
    argangmnt.rIJK.data[0] = (posvel.rIJK.data[1]*posvel.vk-posvel.vj*posvel.rIJK.data[2]);
    argangmnt.rIJK.data[1] = (posvel.rIJK.data[2]*posvel.vi-posvel.rIJK.data[0]*posvel.vk);
    argangmnt.rIJK.data[2] = (posvel.rIJK.data[0]*posvel.vj-posvel.rIJK.data[1]*posvel.vi);
    double argangmntnorm=sqrt(pow2(argangmnt.rIJK.data[0])+pow2(argangmnt.rIJK.data[1])+pow2(argangmnt.rIJK.data[2]));
    koeret.rIJK.data[0]nc = acos(argangmnt.rIJK.data[2]/argangmntnorm);
    nodalvect.rIJK.data[0]=-argangmnt.rIJK.data[1];
    nodalvect.rIJK.data[1]=argangmnt.rIJK.data[0];
    double nodalvectnorm =  atan2(nodalvect.rIJK.data[0],nodalvect.rIJK.data[1]);
    koeret.asc = acos(nodalvect.rIJK.data[0]/nodalvectnorm);
    //if Ny > 0 then 0<asc<180
    //if Ny < 0 then 180<asc<360
    if (nodalvect.rIJK.data[1]<0){
        koeret.asc=2*M_PI-koeret.asc;
    };
    eccvect.p = (1/planet.sgp)*pow2(orbitalvelocity);
    eccvect.q = (-1/planet.sgp)*(posvel.rIJK.data[0]*posvel.vi+posvel.rIJK.data[1]*posvel.rIJK.data[1]+posvel.rIJK.data[2]*posvel.vk);
    eccvect.w = 0;
    eccvect.vp = 0; eccvect.vq = 0; eccvect.vw;
    eci ecceci = PCStoECI(koeret, eccvect);// this may cause an error later
    koeret.ecc = sqrt(pow2(eccvect.p)+pow2(eccvect.q));
    koeret.sma=(pow2(argangmntnorm)/planet.sgp)/(1-pow2(koeret.ecc));
    koeret.aop=acos((nodalvect.rIJK.data[0]*eccvect.p+nodalvect.rIJK.data[1]*eccvect.q)/(nodalvectnorm*koeret.ecc));
    //if ecc_z > 0 then 0<aop<180
    //if ecc_z < 0 then 180<aop<360
    if(ecceci.rIJK.data[2]<0){
        koeret.aop=2*M_PI-koeret.aop;
    };
    koeret.truAnom=acos((koeret.sma*(1-pow2(koeret.ecc))-orbitalradius)/(koeret.ecc*orbitalradius));
    if(posvel.rIJK.data[0]*posvel.vi+posvel.rIJK.data[1]*posvel.vj*posvel.rIJK.data[2]*posvel.vk<0){//dot product of position and velocity
        koeret.truAnom=2*M_PI-koeret.truAnom;
    };
    koeret.eccAnom=atan2(sqrt(1-koeret.ecc)*sin(koeret.truAnom/2),sqrt(1+koeret.ecc)*cos(koeret.truAnom/2));
    koeret.meanAnom=koeret.eccAnom-koeret.ecc*sin(koeret.eccAnom);
    return koeret;
};

ecef ECItoECEF(eci arg,double siderealangle){
    ecef ecefret;
    ecefret.rXYZ.data[0] = cos(siderealangle)*arg.rIJK.data[0]-sin(siderealangle)*arg.rIJK.data[1];
    ecefret.rXYZ.data[1] = sin(siderealangle)*arg.rIJK.data[0]+cos(siderealangle)*arg.rIJK.data[1];
    ecefret.rXYZ.data[2] = arg.rIJK.data[2]; //ignoring precession and nutation
    return ecefret;
};

eci ECEFtoECI(ecef arg, double siderealangle){
    eci eciret;

    eciret.rIJK.data[0]= cos(siderealangle)*arg.rXYZ.data[0]+sin(siderealangle)*arg.rXYZ.data[1];
    eciret.rIJK.data[1]=-sin(siderealangle)*arg.rXYZ.data[0]+cos(siderealangle)*arg.rXYZ.data[1];
    eciret.rIJK.data[2]=arg.rXYZ.data[2]; //ignoring precession and nutation
    return eciret;
};

//the following two are not exactly true axis transforms, and are only used for maneuvers and pointing command handling

eci VNBtoECI(eci posvel, vnb VNBdV){
    eci ECIdV;
    //precalculating as this needs to be run quite fast
    eci rxv, vxrxv;
    rxv.rIJK.data[0]=posvel.rIJK.data[1]*posvel.vk-posvel.rIJK.data[2]*posvel.vj;
    rxv.rIJK.data[1]=posvel.rIJK.data[2]*posvel.vi-posvel.rIJK.data[0]*posvel.vk;
    rxv.rIJK.data[2]=posvel.rIJK.data[0]*posvel.vj-posvel.rIJK.data[1]*posvel.vi;

    vxrxv.rIJK.data[0]=posvel.vj*rxv.rIJK.data[2]-posvel.vk*rxv.rIJK.data[1];
    vxrxv.rIJK.data[1]=posvel.vk*rxv.rIJK.data[0]-posvel.vi*rxv.rIJK.data[2];
    vxrxv.rIJK.data[2]=posvel.vi*rxv.rIJK.data[1]-posvel.vj*rxv.rIJK.data[0];

    double absv, absrxv, absvxrxv;
    absv= sqrt(pow2(posvel.vi)+pow2(posvel.vj)+pow2(posvel.vk));
    absrxv=sqrt(pow2(rxv.rIJK.data[0])+pow2(rxv.rIJK.data[1])+pow2(rxv.rIJK.data[2]));
    absvxrxv=sqrt(pow2(vxrxv.rIJK.data[0])+pow2(vxrxv.rIJK.data[1])+pow2(vxrxv.rIJK.data[2]));



    ECIdV.vi       =(posvel.vi/absv)*VNBdV.v
                    +(rxv.rIJK.data[0]/absrxv)*VNBdV.n
                    +(vxrxv.rIJK.data[0]/absvxrxv)*VNBdV.b;

    ECIdV.vi       =(posvel.vj/absv)*VNBdV.v
                    +(rxv.rIJK.data[1]/absrxv)*VNBdV.n
                    +(vxrxv.rIJK.data[1]/absvxrxv)*VNBdV.b;

    ECIdV.vi       =(posvel.vk/absv)*VNBdV.v
                    +(rxv.rIJK.data[2]/absrxv)*VNBdV.n
                    +(vxrxv.rIJK.data[2]/absvxrxv)*VNBdV.b;
    return ECIdV;
};

vnb ECItoVNB(eci posvel, eci ECIdV){
    vnb VNBdV;
    //precalculating as this needs to be run quite fast
    eci rxv, vxrxv;
    rxv.rIJK.data[0]=posvel.rIJK.data[1]*posvel.vk-posvel.rIJK.data[2]*posvel.vj;
    rxv.rIJK.data[1]=posvel.rIJK.data[2]*posvel.vi-posvel.rIJK.data[0]*posvel.vk;
    rxv.rIJK.data[2]=posvel.rIJK.data[0]*posvel.vj-posvel.rIJK.data[1]*posvel.vi;

    vxrxv.rIJK.data[0]=posvel.vj*rxv.rIJK.data[2]-posvel.vk*rxv.rIJK.data[1];
    vxrxv.rIJK.data[1]=posvel.vk*rxv.rIJK.data[0]-posvel.vi*rxv.rIJK.data[2];
    vxrxv.rIJK.data[2]=posvel.vi*rxv.rIJK.data[1]-posvel.vj*rxv.rIJK.data[0];

    double absv, absrxv, absvxrxv;
    absv= sqrt(pow2(posvel.vi)+pow2(posvel.vj)+pow2(posvel.vk));
    absrxv=sqrt(pow2(rxv.rIJK.data[0])+pow2(rxv.rIJK.data[1])+pow2(rxv.rIJK.data[2]));
    absvxrxv=sqrt(pow2(vxrxv.rIJK.data[0])+pow2(vxrxv.rIJK.data[1])+pow2(vxrxv.rIJK.data[2]));

    VNBdV.v       =(posvel.vi/absv)*ECIdV.vi+(posvel.vj/absv)*ECIdV.vj+(posvel.vk/absv)*ECIdV.vk;

    VNBdV.n       =(rxv.rIJK.data[0]/absrxv)*ECIdV.vi+(rxv.rIJK.data[1]/absrxv)*ECIdV.vj+(rxv.rIJK.data[2]/absrxv)*ECIdV.vk;

    VNBdV.b       =(vxrxv.rIJK.data[0]/absvxrxv)*ECIdV.vi+(vxrxv.rIJK.data[1]/absvxrxv)*ECIdV.vj+(vxrxv.rIJK.data[2]/absvxrxv)*ECIdV.vk;
    return VNBdV;
};

};

    