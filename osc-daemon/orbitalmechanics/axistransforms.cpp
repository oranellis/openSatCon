#include <iostream>
#include <vector>
#include <math.h>
#include <chrono>

#include "../osctypes.hpp"
#include "planet.cpp"

namespace osc{

orbParam meanToTrue(orbParam KOE) {
    // converts mean anomaly to true anomaly
    if (KOE.ecc < 0.2) { //this is a handy calculation to save time for small eccentricities
        KOE.truAnom = KOE.meanAnom + 2 * KOE.ecc * sin(KOE.meanAnom) 
                      +1.25 * pow2(KOE.ecc) * sin(2 * KOE.meanAnom) 
                      - pow3(KOE.ecc) * (0.25 * sin(KOE.meanAnom) - (13/12) * sin(3 * KOE.meanAnom));
    }
    
    else if (KOE.ecc < 1.0) {// newton raphson's method must be used for higher eccentricities, e>1 is a parabolic orbit
        KOE.eccAnom = KOE.meanAnom + ((KOE.ecc * sin(KOE.meanAnom) / (cos(KOE.ecc) - (M_PI_2 - KOE.ecc) * sin(KOE.ecc) + KOE.meanAnom * sin(KOE.ecc))));
        double dE = KOE.eccAnom;
        while(abs(dE) > 10e-10) {
            dE = (KOE.eccAnom - KOE.ecc * sin(KOE.eccAnom) - KOE.meanAnom) / (1 - KOE.ecc * cos(KOE.eccAnom));
            KOE.eccAnom = KOE.eccAnom - dE;
        };
        KOE.truAnom = atan2(sqrt(1 - pow2(KOE.ecc) * sin(KOE.ecc)), cos(KOE.eccAnom) - KOE.ecc);
    };
        return KOE;
};


double greenwichSiderealAngle() { //not working
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
    double julianDays = seconds/1314900;
    gsraret = remainder(2 * M_PI * (0.77905722732640 + 1.00273781191135448 * (julianDays - 2451545)), 2*M_PI);

    return gsraret;
    //this could be made down to sub-second accuracy if possible
};


ecef LLAtoECEF(lla posLLA) {//working
    ecef posECEF;
    // converts the angular position of the satellite to the ECEF co-ordinate system
    // this gives the correct altitude of a non-spherical earth, note that altitude here is above the ground
    double normalDistance = planet.sMa / sqrt(1 - pow2(planet.ecc * sin(posLLA.lat)));

        posECEF.rXYZ.data[0] = (normalDistance + posLLA.alt) * cos(posLLA.lon) * cos(posLLA.lat);
        posECEF.rXYZ.data[1] = (normalDistance + posLLA.alt) * sin(posLLA.lon) * cos(posLLA.lat);
        posECEF.rXYZ.data[2] = (normalDistance * (1 - pow2(planet.ecc)) + posLLA.alt) * sin(posLLA.lat);

    return posECEF;
};

lla ECEFtoLLA(ecef posECEF) {//working
    //returns ground position of sub satellite point and satellite altitude from ECEF co-ords
    lla llaret;
    double secondeccentricity = planet.ecc / sqrt(1 - pow2(planet.ecc));
    double p = sqrt(pow2(posECEF.rXYZ.data[0]) + pow2(posECEF.rXYZ.data[1]));
    double theta = atan2(posECEF.rXYZ.data[2] * planet.sMa, p * planet.sma);

        llaret.lon = atan2(posECEF.rXYZ.data[1], posECEF.rXYZ.data[0]);
        llaret.lat = atan2(posECEF.rXYZ.data[2] + pow2(secondeccentricity) * planet.sma * pow3(sin(theta)),
                     p - pow2(planet.ecc) * planet.sMa * pow3(cos(theta)));

    double normaldistance = planet.sMa / sqrt(1 - pow2(planet.ecc*sin(llaret.lat)));

        llaret.alt = (p / cos(llaret.lat)) - normaldistance;

    return llaret;
};

ned ECEFtoNED(ecef posECEF, lla refLLA) { //working
    // returns vector from satellite to ground station in North East Down reference frame
    ned posNED;
    ecef refECEF = LLAtoECEF(refLLA);
    ecef relECEF;
    relECEF.rXYZ = posECEF.rXYZ.operator-(refECEF.rXYZ);

        posNED.rNED.data[0] = -sin(refLLA.lat) * cos(refLLA.lon) * relECEF.rXYZ.data[0]
                              -sin(refLLA.lat) * sin(refLLA.lon) * relECEF.rXYZ.data[1]
                              +cos(refLLA.lat)                   * relECEF.rXYZ.data[2];

        posNED.rNED.data[1] = -sin(refLLA.lon)                   * relECEF.rXYZ.data[0]
                              +cos(refLLA.lon)                   * relECEF.rXYZ.data[1]
                              +0;

        posNED.rNED.data[2] = -cos(refLLA.lat) * cos(refLLA.lon) * relECEF.rXYZ.data[0]
                              -cos(refLLA.lat) * sin(refLLA.lon) * relECEF.rXYZ.data[1]
                              -sin(refLLA.lat)                   * relECEF.rXYZ.data[2];

    return posNED;
};

ecef NEDtoECEF(ned posNED, lla refLLA) {//working
    // returns ECEF position from NED vector
    ecef relECEF;
    ecef refECEF = LLAtoECEF(refLLA);
    ecef posECEF;

        relECEF.rXYZ.data[0] = -sin(refLLA.lat) * cos(refLLA.lon) * posNED.rNED.data[0]
                               -sin(refLLA.lon)                   * posNED.rNED.data[1]
                               -cos(refLLA.lat) * cos(refLLA.lon) * posNED.rNED.data[2];

        relECEF.rXYZ.data[1] = -sin(refLLA.lat) * sin(refLLA.lon) * posNED.rNED.data[0]
                               +cos(refLLA.lon)                   * posNED.rNED.data[1]
                               -cos(refLLA.lat) * sin(refLLA.lon) * posNED.rNED.data[2];

        relECEF.rXYZ.data[2] = cos(refLLA.lat)                    * posNED.rNED.data[0]
                               +0
                               -sin(refLLA.lat)                   * posNED.rNED.data[2];

    posECEF.rXYZ = relECEF.rXYZ.operator+(refECEF.rXYZ);
    return posECEF;   
};

ecef EARtoECEF(ear posEAR, lla refLLA) {
    // returns ECEF position from ground station Elevation Azimuth Range tracking
    // this method goes through SEZ co-ordinates, but seems to be the best method
    ecef posECEF;
    ecef relECEF;
    ecef refECEF = LLAtoECEF(refLLA);

        relECEF.rXYZ.data[0] = sin(refLLA.lat) * cos(refLLA.lon)  * (-posEAR.r * cos(posEAR.e) * cos(posEAR.a))
                               -sin(refLLA.lon)                   * (posEAR.r * cos(posEAR.e) * sin(posEAR.a))
                               +cos(refLLA.lat) * cos(refLLA.lon) * (posEAR.r * sin(posEAR.e));
        
        relECEF.rXYZ.data[1] = sin(refLLA.lat) * sin(refLLA.lon)  * (-posEAR.r * cos(posEAR.e) * cos(posEAR.a))
                               +cos(refLLA.lon)                   * (posEAR.r * cos(posEAR.e) * sin(posEAR.a))
                               +cos(refLLA.lat) * sin(refLLA.lon) * (posEAR.r * sin(posEAR.e));

        relECEF.rXYZ.data[2] = -cos(refLLA.lat)                   * (-posEAR.r * cos(posEAR.e) * cos(posEAR.a))
                               +0
                               +sin(refLLA.lat)                   * (posEAR.r * sin(posEAR.e));

    posECEF.rXYZ=relECEF.rXYZ.operator+(refECEF.rXYZ);

    return posECEF;
};

ecef ENUtoECEF(enu posENU, lla refLLA) {
    // returns the ECEF position of a satellite tracked by a ground station in East North Up
    ecef posECEF;
    ecef relECEF;
    ecef refECEF = LLAtoECEF(refLLA);

        relECEF.rXYZ.data[0] = -sin(refLLA.lon)                   * posENU.rENU.data[0]
                               -sin(refLLA.lat) * cos(refLLA.lon) * posENU.rENU.data[1]
                               +cos(refLLA.lat) * cos(refLLA.lon) * posENU.rENU.data[2];

        relECEF.rXYZ.data[1] = cos(refLLA.lon)                    * posENU.rENU.data[0]
                               -sin(refLLA.lat) * sin(refLLA.lon) * posENU.rENU.data[1]
                               +cos(refLLA.lat) * sin(refLLA.lon) * posENU.rENU.data[2];

        relECEF.rXYZ.data[2] = 0
                               +cos(refLLA.lat)                   * posENU.rENU.data[1]
                               +sin(refLLA.lat)                   * posENU.rENU.data[2];

    posECEF.rXYZ.operator+(refECEF.rXYZ);

    return posECEF;
};

enu ECEFtoENU(ecef posECEF, lla refLLA){
    // returns the ENU values seen by a ground station of a satellite at the given ECEF position
    enu posENU;
    ecef refECEF = LLAtoECEF(refLLA);
    ecef relECEF;
    relECEF.rXYZ = posECEF.rXYZ.operator-(refECEF.rXYZ);

        posENU.rENU.data[1] = -sin(refLLA.lon)                   * (relECEF.rXYZ.data[0])
                               +cos(refLLA.lon)                  * (relECEF.rXYZ.data[1])
                               +0;

        posENU.rENU.data[1] = -sin(refLLA.lat) * cos(refLLA.lon) * relECEF.rXYZ.data[0]
                              -sin(refLLA.lat) * sin(refLLA.lon) * relECEF.rXYZ.data[1]
                              +cos(refLLA.lat)                   * relECEF.rXYZ.data[2];

        posENU.rENU.data[1] = cos(refLLA.lat) * cos(refLLA.lon)  * relECEF.rXYZ.data[0]
                              +cos(refLLA.lat) * sin(refLLA.lon) * relECEF.rXYZ.data[1]
                              +sin(refLLA.lat)                   * relECEF.rXYZ.data[2];

    return posENU;
};

ear ECEFtoEAR(ecef posECEF, lla refLLA) {
    // returns the EAR values seen by a ground station of a satellite at the given ECEF position
    ear posEAR;
    ecef refECEF = LLAtoECEF(refLLA);
    ecef relECEF;
    relECEF.rXYZ = posECEF.rXYZ.operator-(refECEF.rXYZ);
    enu relENU = ECEFtoENU(relECEF, refLLA);


        posEAR.r = relENU.rENU.mag();
        posEAR.e = asin(relENU.rENU.data[2] / posEAR.r);
        posEAR.a = atan2(relENU.rENU.data[0], relENU.rENU.data[1]);
    
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

ecef SEZtoECEF(thcs posSEZ, lla refLLA) {
    // returns the ECEF position of a satellite tracked by a ground station in the Topocentric Horizon Co-ordinate System
    ecef posECEF;
    ecef relECEF;
    ecef refECEF = LLAtoECEF(refLLA);
    refLLA.lat = M_PI_2 - refLLA.lat; // rotates to horizon

        relECEF.rXYZ.data[0] = cos(refLLA.lat) * cos(refLLA.lon)  * posSEZ.rSEZ.data[0]
                               -sin(refLLA.lon)                   * posSEZ.rSEZ.data[1]
                               +sin(refLLA.lat) * cos(refLLA.lat) * posSEZ.rSEZ.data[2];

        relECEF.rXYZ.data[1] = cos(refLLA.lat) * sin(refLLA.lon)  * posSEZ.rSEZ.data[0]
                               +cos(refLLA.lon)                   * posSEZ.rSEZ.data[1]
                               +sin(refLLA.lat) * cos(refLLA.lon) * posSEZ.rSEZ.data[2];

        relECEF.rXYZ.data[2] = -sin(refLLA.lat)                   * posSEZ.rSEZ.data[0]
                               +0
                               +cos(refLLA.lat)                   * posSEZ.rSEZ.data[2];

    posECEF.rXYZ.operator+(refECEF.rXYZ);

    return posECEF;
};

thcs ECEFtoSEZ(ecef posECEF, lla refLLA) {
    // returns the SEZ values seen by a ground station of a satellite at the given ECEF position
    thcs posSEZ;
    ecef refECEF = LLAtoECEF(refLLA);
    ecef relECEF;
    relECEF.rXYZ = posECEF.rXYZ.operator-(refECEF.rXYZ);
    refLLA.lat = M_PI_2 - refLLA.lat; // rotates to horizon

        posSEZ.rSEZ.data[0] = cos(refLLA.lat) * cos(refLLA.lon)  * relECEF.rXYZ.data[0]
                              +cos(refLLA.lat) * sin(refLLA.lon) * relECEF.rXYZ.data[1]
                              -sin(refLLA.lat)                   * relECEF.rXYZ.data[2];

        posSEZ.rSEZ.data[1] = -sin(refLLA.lon)                   * relECEF.rXYZ.data[0]
                              +cos(refLLA.lon)                   * relECEF.rXYZ.data[1]
                              +0;

        posSEZ.rSEZ.data[2] = sin(refLLA.lat) * cos(refLLA.lon)  * relECEF.rXYZ.data[0]
                              +sin(refLLA.lat) * sin(refLLA.lon) * relECEF.rXYZ.data[1]
                              +cos(refLLA.lat)                   * relECEF.rXYZ.data[2];

    return posSEZ;
};

eci PCStoECI(orbParam KOE, pcs posvelPCS) { //working
    eci posvelECI;

        posvelECI.rIJK.data[0] = (cos(KOE.aop) * cos(KOE.asc) - sin(KOE.asc) * sin(KOE.aop) * cos(KOE.inc))   * posvelPCS.rPCS.data[0]
                                 +(-cos(KOE.asc) * sin(KOE.aop) - sin(KOE.asc) * cos(KOE.aop) * cos(KOE.inc)) * posvelPCS.rPCS.data[1]
                                 +(sin(KOE.asc) * sin(KOE.inc))                                               * posvelPCS.rPCS.data[2];

        posvelECI.rIJK.data[1] = (sin(KOE.asc) * cos(KOE.aop) + cos(KOE.asc) * sin(KOE.aop) * cos(KOE.inc))   * posvelPCS.rPCS.data[0]
                                 +(-sin(KOE.asc) * sin(KOE.aop) + cos(KOE.asc) * cos(KOE.aop) * cos(KOE.inc)) * posvelPCS.rPCS.data[1]
                                 +(-cos(KOE.asc) * sin(KOE.inc))                                              * posvelPCS.rPCS.data[2];

        posvelECI.rIJK.data[2] = sin(KOE.aop) * sin(KOE.inc)                                                  * posvelPCS.rPCS.data[0]
                                 +cos(KOE.aop) * sin(KOE.inc)                                                 * posvelPCS.rPCS.data[1]
                                 +cos(KOE.inc)                                                                * posvelPCS.rPCS.data[2];



        posvelECI.vIJK.data[0] = (cos(KOE.aop) * cos(KOE.asc) - sin(KOE.asc) * sin(KOE.aop) * cos(KOE.inc))   * posvelPCS.vPCS.data[0]
                                 +(-cos(KOE.asc) * sin(KOE.aop) - sin(KOE.asc) * cos(KOE.aop) * cos(KOE.inc)) * posvelPCS.vPCS.data[1]
                                 +(sin(KOE.asc) * sin(KOE.inc))                                               * posvelPCS.vPCS.data[2];

        posvelECI.vIJK.data[1] = (sin(KOE.asc) * cos(KOE.aop) + cos(KOE.asc) * sin(KOE.aop) * cos(KOE.inc))   * posvelPCS.vPCS.data[0]
                                 +(-sin(KOE.asc) * sin(KOE.aop) + cos(KOE.asc) * cos(KOE.aop) * cos(KOE.inc)) * posvelPCS.vPCS.data[1]
                                 +(-cos(KOE.asc) * sin(KOE.inc))                                              * posvelPCS.vPCS.data[2];

        posvelECI.vIJK.data[2] = sin(KOE.aop) * sin(KOE.inc)                                                  * posvelPCS.vPCS.data[0]
                                 +cos(KOE.aop) * sin(KOE.inc)                                                 * posvelPCS.vPCS.data[1]
                                 +cos(KOE.inc)                                                                * posvelPCS.vPCS.data[2];

    return posvelECI;
};

pcs KOEtoPCS(orbParam KOE) { //working
    // gives the PCS position and velocity of a satellite given its Kepler Orbital Elements

    pcs posvelPCS;
    double orbitRadius = (KOE.sma * (1 - pow2(KOE.ecc)))/(1 + KOE.ecc * cos(KOE.truAnom));

        posvelPCS.rPCS.data[0] = orbitRadius * cos(KOE.truAnom); 
        posvelPCS.rPCS.data[1] = orbitRadius * sin(KOE.truAnom);
        posvelPCS.rPCS.data[2] = 0;

    double p        = KOE.sma * (1 - pow2(KOE.ecc));
    double thetaDot = sqrt(planet.sgp * p) / pow2(orbitRadius);
    double rDot     = sqrt(planet.sgp / p) * KOE.ecc * sin(KOE.truAnom);

        posvelPCS.vPCS.data[0] = rDot * cos(KOE.truAnom) - orbitRadius * sin(KOE.truAnom) * thetaDot;
        posvelPCS.vPCS.data[1] = rDot * sin(KOE.truAnom) + orbitRadius * cos(KOE.truAnom) * thetaDot;
        posvelPCS.vPCS.data[2] = 0;

    return posvelPCS;
};

orbParam ECItoKOE(eci posvelECI) {
    orbParam KOE;
    eci h;                           // angular momentum vector
    eci N;                           // ascending node vector
    eci e;                           // eccentricity vector
    double r = posvelECI.rIJK.mag(); // orbital radius
    double v = posvelECI.vIJK.mag(); // orbital velocity

    h.rIJK = posvelECI.rIJK.cross(posvelECI.vIJK);

    double normH=h.rIJK.mag();

    KOE.inc = acos(h.rIJK.data[2] / normH);

    N.rIJK.data[0] = -h.rIJK.data[1];
    N.rIJK.data[1] = h.rIJK.data[0];

    double normN = N.rIJK.mag();
    KOE.asc = acos(N.rIJK.data[0] / normN);
    //if Ny > 0 then 0<asc<180
    //if Ny < 0 then 180<asc<360
    if (N.rIJK.data[1] < 0) {
        KOE.asc = 2 * M_PI - KOE.asc;
    };

    double vSquaredMinusMuOverR = (pow2(v) - (planet.sgp / r));
    eci vSquaredMinusMuOverRtimesR, rDotVtimesV, vSquaredMinusMuOverRtimesRminusrDotVtimesV; //fix this horror
    double rDotV;

    vSquaredMinusMuOverRtimesR.rIJK                 = posvelECI.rIJK.operator*(vSquaredMinusMuOverR);
    rDotV                                           = posvelECI.rIJK.dot(posvelECI.vIJK);
    rDotVtimesV.rIJK                                = posvelECI.vIJK.operator*(rDotV);
    vSquaredMinusMuOverRtimesRminusrDotVtimesV.rIJK = vSquaredMinusMuOverRtimesR.rIJK.operator-(rDotVtimesV.rIJK);

    e.rIJK                                          = vSquaredMinusMuOverRtimesRminusrDotVtimesV.rIJK.operator/(planet.sgp);

    KOE.ecc = e.rIJK.mag();
    KOE.sma = (pow2(normH) / planet.sgp) / (1 - pow2(KOE.ecc));
    KOE.aop = acos((N.rIJK.dot(e.rIJK) / (normN * KOE.ecc)));

    //if ecc_z > 0 then 0<aop<180
    //if ecc_z < 0 then 180<aop<360

    if(e.rIJK.data[2] < 0) {
        KOE.aop = 2 * M_PI - KOE.aop;
    };

    KOE.truAnom = acos((e.rIJK.dot(posvelECI.rIJK) / (r*KOE.ecc)));

    if(posvelECI.rIJK.dot(posvelECI.vIJK) < 0) {//dot product of position and velocity
        KOE.truAnom = 2 * M_PI - KOE.truAnom;
    };

    KOE.eccAnom  = atan2(sqrt(1 - KOE.ecc) * sin(KOE.truAnom / 2), sqrt(1 + KOE.ecc) * cos(KOE.truAnom / 2));
    KOE.meanAnom = KOE.eccAnom - KOE.ecc * sin(KOE.eccAnom);

    return KOE;
};

ecef ECItoECEF(eci posvelECI,double siderealAngle) {//working
    ecef posvelECEF;

        posvelECEF.rXYZ.data[0] = cos(siderealAngle) * posvelECI.rIJK.data[0] + sin(siderealAngle) * posvelECI.rIJK.data[1];
        posvelECEF.rXYZ.data[1] = -sin(siderealAngle) * posvelECI.rIJK.data[0] + cos(siderealAngle) * posvelECI.rIJK.data[1];
        posvelECEF.rXYZ.data[2] =                                                                    posvelECI.rIJK.data[2]; //ignoring precession and nutation

    return posvelECEF;
};

eci ECEFtoECI(ecef posvelECEF, double siderealAngle) {//working
    eci posvelECI;

        posvelECI.rIJK.data[0] = cos(-siderealAngle) * posvelECEF.rXYZ.data[0] + sin(-siderealAngle)  * posvelECEF.rXYZ.data[1];
        posvelECI.rIJK.data[1] = -sin(-siderealAngle) * posvelECEF.rXYZ.data[0] + cos(-siderealAngle) * posvelECEF.rXYZ.data[1];
        posvelECI.rIJK.data[2] =                                                                        posvelECEF.rXYZ.data[2]; //ignoring precession and nutation
    
    return posvelECI;
};

//the following two are not exactly true axis transforms, and are only used for maneuvers and pointing command handling

eci VNBtoECI(eci posvelECI, vnb VNBdV) {
    eci ECIdV;
    //precalculating as this needs to be run quite fast
    eci rxv, vxrxv;
    rxv.rIJK   = posvelECI.rIJK.cross(posvelECI.vIJK);
    vxrxv.rIJK = posvelECI.vIJK.cross(rxv.rIJK);

    double absv, absrxv, absvxrxv;
    absv     = posvelECI.vIJK.mag();
    absrxv   = rxv.rIJK.mag();
    absvxrxv = vxrxv.rIJK.mag();

        ECIdV.vIJK.data[0] = (posvelECI.vIJK.data[0] / absv)  * VNBdV.vVNB.data[0]
                             +(rxv.rIJK.data[0] / absrxv)     * VNBdV.vVNB.data[1]
                             +(vxrxv.rIJK.data[0] / absvxrxv) * VNBdV.vVNB.data[2];

        ECIdV.vIJK.data[1] = (posvelECI.vIJK.data[1] / absv)  * VNBdV.vVNB.data[0]
                             +(rxv.rIJK.data[1] / absrxv)     * VNBdV.vVNB.data[1]
                             +(vxrxv.rIJK.data[1] / absvxrxv) * VNBdV.vVNB.data[2];

        ECIdV.vIJK.data[2] = (posvelECI.vIJK.data[2] / absv)  * VNBdV.vVNB.data[0]
                             +(rxv.rIJK.data[2] / absrxv)     * VNBdV.vVNB.data[1]
                             +(vxrxv.rIJK.data[2] / absvxrxv) * VNBdV.vVNB.data[2];

    return ECIdV;
};

vnb ECItoVNB(eci posvelECI, eci ECIdV) {
    vnb VNBdV;
    //precalculating as this needs to be run quite fast
    eci rxv, vxrxv;
    rxv.rIJK   = posvelECI.rIJK.cross(posvelECI.vIJK);
    vxrxv.rIJK = posvelECI.vIJK.cross(rxv.rIJK);

    double absv, absrxv, absvxrxv;
    absv     = posvelECI.vIJK.mag();
    absrxv   = rxv.rIJK.mag();
    absvxrxv = vxrxv.rIJK.mag();

        VNBdV.vVNB.data[0] = (posvelECI.vIJK.data[0] / absv)  * ECIdV.vIJK.data[0]
                             +(posvelECI.vIJK.data[1] / absv) * ECIdV.vIJK.data[1]
                             +(posvelECI.vIJK.data[2] / absv) * ECIdV.vIJK.data[2];
        
        VNBdV.vVNB.data[1] = (rxv.rIJK.data[0] / absrxv)      * ECIdV.vIJK.data[0]
                             +(rxv.rIJK.data[1] / absrxv)     * ECIdV.vIJK.data[1]
                             +(rxv.rIJK.data[2] / absrxv)     * ECIdV.vIJK.data[2];

        VNBdV.vVNB.data[2] = (vxrxv.rIJK.data[0] / absvxrxv)  * ECIdV.vIJK.data[0]
                             +(vxrxv.rIJK.data[1] /absvxrxv)  * ECIdV.vIJK.data[1]
                             +(vxrxv.rIJK.data[2] /absvxrxv)  * ECIdV.vIJK.data[2];
                             
    return VNBdV;
};

};

    