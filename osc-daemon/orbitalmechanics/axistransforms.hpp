#ifndef AXISTRANSFORMS_H
#define AXISTRANSFORMS_H

#include <iostream>
#include <vector>
#include <math.h>
#include <chrono>

#include "../osctypes.hpp"
#include "planet.hpp"

namespace osc{

    /** greenwichSiderealAngle() calculates the angle required for the transformation between ECI and ECEF
    accounts for the rotation of the Earth using the number of Julian days
    since the J2000 epoch, using an equation from the IAU
    */
inline double greenwichSiderealAngle() { //working

    double gsraret;
    time_t timer;
    //for Oran
    struct tm J2000 = {0};
    double seconds;

    J2000.tm_hour = 11;  J2000.tm_min = 58; J2000.tm_sec = 55;
    J2000.tm_year = 100; J2000.tm_mon = 0; J2000.tm_mday = 1;

    time(&timer);

    seconds = difftime(timer,mktime(&J2000));
    double julianDays = seconds/86400;
    
    gsraret = fmod(2 * M_PI * (0.77905722732640 + 1.00273781191135448 * (julianDays)), 2 * M_PI);

    return gsraret;
    //this could be made down to sub-second accuracy if possible
};

    /** LLAtoECEF(posLLA) converts the latitude, longitude, and altitude of an object to the ECEF co-ordinate system
    this uses the correct altitude of a non-spherical earth, note that altitude here is above the ground
    @param[in] posLLA position in LLA coordinates
    */
inline ecef LLAtoECEF(lla posLLA) {
    ecef posECEF;
    double normalDistance = planet.sMa / sqrt(1 - pow2(planet.ecc * sin(posLLA.lat)));

    posECEF.rXYZ.data[0] = (normalDistance + posLLA.alt) * cos(posLLA.lon) * cos(posLLA.lat);
    posECEF.rXYZ.data[1] = (normalDistance + posLLA.alt) * sin(posLLA.lon) * cos(posLLA.lat);
    posECEF.rXYZ.data[2] = (normalDistance * (1 - pow2(planet.ecc)) + posLLA.alt) * sin(posLLA.lat);

    return posECEF;
};

    /** \ECEFtoLLA(posECEF) returns ground position of sub satellite point and satellite altitude from ECEF co-ords
    @param[in] posECEF input ECEF coordinate position
    */

inline lla ECEFtoLLA(ecef posECEF) {
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

    /** ECEFtoNED(posECEF, refLLA) returns a vector from satellite to ground station in North East Down reference frame
    this axis frame intuitively forms a local tangent plane, and has good axes for 
    visualisation, unlike ECI or ECEF
    @param[in] posECEF ECEF coordinate position
    @param[in] refLLA reference LLA position

    */
inline ned ECEFtoNED(ecef posECEF, lla refLLA) {
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

    /** NEDtoECEF(posNED, refLLA) returns an ECEF from an input NED.
    Inverse of ECEFtoNED, can be used for satellite measuring against a ground object
    @param[in] posNED NED coordinate position
    @param[in] refLLA reference LLA

    */
inline ecef NEDtoECEF(ned posNED, lla refLLA) {
    // 
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

    /** EARtoECEF(posEAR, refLLA) returns ECEF position from ground station Elevation Azimuth Range tracking
    this method goes through SEZ co-ordinates, but seems to be the best method
    @param[in] posEAR position EAR coordinates
    @param[in] refLLA reference LLA

    */
inline ecef EARtoECEF(ear posEAR, lla refLLA) {
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

    /** ENUtoECEF(posENU, refLLA) returns the ECEF position of a satellite tracked by a ground station in East North Up
    this is simply an axis shift of NED, but is used when ground stations track satellites
    instead of vice versa
    @param[in] posENU NED coordinate position
    @param[in] refLLA reference LLA

    */
inline ecef ENUtoECEF(enu posENU, lla refLLA) {

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

    /** ECEFtoENU(posECEF, refLLA) returns the ENU values seen by a ground station of a satellite at the given ECEF position
    used as an intermediate step for the following transform
    @param[in] posECEF ECEF coordinate position
    @param[in] refLLA reference LLA

    */
inline enu ECEFtoENU(ecef posECEF, lla refLLA){
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

    /** ECEFtoEAR(posECEF, refLLA) returns the EAR values seen by a ground station of a satellite at the given ECEF position
    @param[in] posECEF ECEF coordinate position
    @param[in] refLLA reference LLA
    */
inline ear ECEFtoEAR(ecef posECEF, lla refLLA) {
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

    /** SEZtoECEF(posSEZ, refLLA) returns the ECEF position of a satellite tracked by a ground station in the Topocentric Horizon Co-ordinate System
    @param[in] posSEZ SEZ coordinate position
    @param[in] refLLA reference LLA
    */
inline ecef SEZtoECEF(thcs posSEZ, lla refLLA) {
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

    /** ECEFtoSEZ(posECEF, refLLA) returns the SEZ values seen by a ground station of a satellite at the given ECEF position
    @param[in] posECEF NED coordinate position
    @param[in] refLLA reference LLA
    */
inline thcs ECEFtoSEZ(ecef posECEF, lla refLLA) {
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

    /** PCStoECI(KOE, posvelPCS) returns an ECI from an input PCS.
    This is the second step for getting ECI positon and velocity from Keplerian elements,
    by using a co-ordinate system with two vectors in the plane of the orbit 
    @param[in] KOE KOE of craft
    @param[in] posvelPCS position and velocity in PCS coordinates
    
    */
inline eci PCStoECI(orbParam KOE, pcs posvelPCS) {
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

    /** KOEtoPCS(KOE) gives the PCS position and velocity of a satellite given its Keplerian Orbital Elements
    used to determine the ECI and ECEF position from Keplerian elements
    @param[in] KOE current KOE

   
    */
inline pcs KOEtoPCS(orbParam KOE) {
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
    /** ECItoKOE(posvelECI) this returns the kepler orbital elements of an orbit from a satellite tracked at a single position and velocity
    which allows orbits to be determined from a single detection of a satellite
    this can also be used tosee how a burn in ECI co-ordinates affects the orbital elements,
    by increasing the velocity vector by the deltaV in the ECI frame
    @param[in] posvelECI ECI coordinate position and velocity

    */
inline orbParam ECItoKOE(eci posvelECI) {

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

    if(e.rIJK.data[2] < 0) {
        KOE.aop = 2 * M_PI - KOE.aop;
    };

    KOE.truAnom = acos((e.rIJK.dot(posvelECI.rIJK) / (r*KOE.ecc)));

    if(posvelECI.rIJK.dot(posvelECI.vIJK) < 0) {//dot product of position and velocity
        KOE.truAnom = 2 * M_PI - KOE.truAnom;
    };

    return KOE;
};

    /** ECItoECEF(posvelECI, siderealAngle) changes from the inertial ECI frame to the non-enertial ECEF frame
    the ECEF frame is more useful for ground positions, whereas the 
    ECI frame is more useful for orbits and celestial objects
    @param[in] posvelECI position and velocity in ECI coordinates
    @param[in] siderealAngle sidereal angle

    */
inline ecef ECItoECEF(eci posvelECI,double siderealAngle) {
    ecef posvelECEF;

        posvelECEF.rXYZ.data[0] = cos(siderealAngle) * posvelECI.rIJK.data[0] + sin(siderealAngle) * posvelECI.rIJK.data[1];
        posvelECEF.rXYZ.data[1] = -sin(siderealAngle) * posvelECI.rIJK.data[0] + cos(siderealAngle) * posvelECI.rIJK.data[1];
        posvelECEF.rXYZ.data[2] =                                                                    posvelECI.rIJK.data[2]; //ignoring precession and nutation

    return posvelECEF;
};

    /** ECEFtoECI(posvelECEF, siderealAngle) rotation in the opposite direction. Complicated effects of precession and nutation have been ignored,
    but can be simply implemented into the transforms
    @param[in] posvelECEF position and velocity in ECEF coordinates
    @param[in] siderealAngle siderealAngle value

    */
inline eci ECEFtoECI(ecef posvelECEF, double siderealAngle) {

    eci posvelECI;

        posvelECI.rIJK.data[0] = cos(-siderealAngle) * posvelECEF.rXYZ.data[0] + sin(-siderealAngle)  * posvelECEF.rXYZ.data[1];
        posvelECI.rIJK.data[1] = -sin(-siderealAngle) * posvelECEF.rXYZ.data[0] + cos(-siderealAngle) * posvelECEF.rXYZ.data[1];
        posvelECI.rIJK.data[2] =                                                                        posvelECEF.rXYZ.data[2]; //ignoring precession and nutation
    
    return posvelECI;
};



    /** VNBtoECI(posvelECEI, VNBdV) is not exactly a true axis transform, and is only used for maneuvers and pointing command handling
    the second called in value in a vector relative to the spacecraft, in the Velocity Normal Bi-Normal frame
    the first called in value is the ECI position and velocity at the time of the transform, and it is required to set up
    the transformation matrix. VNB was chosen as it is not only very intuitive to use (+Velocity to increase orbit, -Velocity 
    to decrease), but it also allows for easy integrated plane changes during another planned orbital maneuver. In comparison
    to other satellite-based reference frames, VNB is orthogonal, allowing the inverse to be easily found
    @param[in] posvelECI position and velocity in ECI coordinates
    @param[in] VNBdV deltaV in VNB coordinate system

    */
inline eci VNBtoECI(eci posvelECI, vnb VNBdV) {
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

    /** ECItoVNB(posvelECI, ECIdV) is not exactly a true axis transform, and is only used for maneuvers and pointing command handling
    the second called in value in a vector relative to the spacecraft, in the Earth Centred Inertial frame
    the first called in value is the ECI position and velocity at the time of the transform, and it is required to set up
    the transformation matrix. VNB was chosen as it is not only very intuitive to use (+Velocity to increase orbit, -Velocity 
    to decrease), but it also allows for easy integrated plane changes during another planned orbital maneuver. In comparison
    to other satellite-based reference frames, VNB is orthogonal, allowing the inverse to be easily found *
    @param[in] posvelECI ECI coordinate position and velocity
    @param[in] ECIdV deltaV in ECI coordinates 
    */
inline vnb ECItoVNB(eci posvelECI, eci ECIdV) {
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
        //these could have been replaced with vector maths but it seemed wise to 
        //keep both matrices similar in form as a visual aid
                             
    return VNBdV;
};

};

#endif