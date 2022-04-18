#include <iostream>
#include <vector>
#include <math.h>
#include <chrono>

#include "osctypes.hpp"
#include "planet.cpp"

namespace osc{

double meantotrue(orbparam KOE){
    // converts mean anomaly to true anomaly
        if (KOE.ecc < 0.2){ //this is a handy calculation to save time for small eccentricities
            KOE.truanom = KOE.meananom+2*KOE.ecc*sin(KOE.meananom)+1.25*pow2(KOE.ecc)*sin(2*KOE.meananom)-pow3(KOE.ecc)*(0.25*sin(KOE.meananom)-(13/12)*sin(3*KOE.meananom));
        } else if(KOE.ecc < 1.0){// newton raphson's method must be used for higher eccentricities, e>1 is a parabolic orbit
            KOE.eccanom=KOE.meananom+((KOE.ecc*sin(KOE.meananom)/(cos(KOE.ecc)-(M_PI_2-KOE.ecc)*sin(KOE.ecc)+KOE.meananom*sin(KOE.ecc))));
            double dE=KOE.eccanom;
            while(abs(dE)>10e-10){
                dE=(KOE.eccanom-KOE.ecc*sin(KOE.eccanom)-KOE.meananom)/(1-KOE.ecc*cos(KOE.eccanom));
                KOE.eccanom=KOE.eccanom-dE;
            };
            KOE.truanom=atan2(sqrt(1-pow2(KOE.ecc)*sin(KOE.ecc)),cos(KOE.eccanom)-KOE.ecc);
        };
};


double greenwichsiderealangle(){
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
    double juliandays=seconds/60/60/365.25;
    gsraret=2*M_PI*(0.77905722732640+1.00273781191135448*(juliandays-2451545));
    return gsraret;
    //this could be made down to sub-second accuracy if possible
};


ecef LLAtoECEF(lla arg) {
    ecef ecefret;
    // converts the angular position of the satellite to the ECEF co-ordinate system
    // this gives the correct altitude of a non-spherical earth, note that altitude here is above the ground
    double normaldistance=planet.sMa/sqrt(1-(planet.ecc*sin(arg.lat)));
    ecefret.x=(normaldistance+arg.alt)*cos(arg.lon)*cos(arg.lat);
    ecefret.y=(normaldistance+arg.alt)*sin(arg.lon)*cos(arg.lat);
    ecefret.z=(normaldistance*(1-pow2(planet.ecc))+arg.alt)*sin(arg.lat);
    return ecefret;
};

lla ECEFtoLLA(ecef arg){
    //returns ground position of sub satellite point and satellite altitude from ECEF co-ords
    lla llaret;
    double secondeccentricity=planet.ecc/sqrt(1-pow2(planet.ecc));
    double p = sqrt(pow2(arg.x)+pow2(arg.y));
    double theta = atan2(arg.z*planet.sMa,p*planet.sma);
        llaret.lon = atan2(arg.y,arg.x);
        llaret.lat = atan2(arg.z+pow2(secondeccentricity)*planet.sma*pow3(sin(theta)),p-pow2(planet.ecc)*planet.sMa*pow3(sin(theta)));
    double normaldistance = planet.sMa/sqrt(1-pow2(planet.ecc)*pow2(sin(llaret.lat)));
        llaret.alt = (p/cos(llaret.lat))-normaldistance;
    return llaret;
};

ned ECEFtoNED(ecef satpos, lla reflla){
    // returns vector from satellite to ground station in North East Down reference frame
    ned nedret;
    ecef refpos = LLAtoECEF(reflla);
    ecef relpos;

        relpos.x=satpos.x-refpos.x;
        relpos.y=satpos.y-refpos.y;
        relpos.z=satpos.z-refpos.z;

    lla reflla = ECEFtoLLA(refpos);

    nedret.n=   -sin(reflla.lat)*cos(reflla.lon)*relpos.x+
                sin(reflla.lat)*sin(reflla.lon)*relpos.y+
                cos(reflla.lat)*relpos.z;

    nedret.e=   -sin(reflla.lon)*relpos.x+
                cos(reflla.lon)*relpos.y+
                0;

    nedret.d=   -cos(reflla.lat)*cos(reflla.lon)*relpos.x+
                -cos(reflla.lat)*sin(reflla.lon)*relpos.y+
                -sin(reflla.lat)*relpos.z;

    return nedret;
};

ecef NEDtoECEF(ned satpos, lla reflla){
    // returns ECEF position from NED vector
    ecef ecefret;

    ecefret.x=  -cos(reflla.lon)*sin(reflla.lat)*satpos.n
                -sin(reflla.lon)*satpos.e
                -cos(reflla.lon)*cos(reflla.lat)*satpos.d;

    ecefret.y=  -sin(reflla.lon)*sin(reflla.lat)*satpos.n
                +cos(reflla.lon)*satpos.e
                -sin(reflla.lon)*cos(reflla.lat)*satpos.d;

    ecefret.z=  cos(reflla.lat)*satpos.n
                +0
                -sin(reflla.lat)*satpos.d;
    return ecefret;   
};

ecef EARtoECEF(ear satpos, lla reflla){
    // returns ECEF position from ground station Elevation Azimuth Range tracking
    // this method goes through SEZ co-ordinates, but seems to be the best method
    ecef ecefret;
    ecef refpos = LLAtoECEF(reflla);
    ecefret.x=  sin(reflla.lat)*cos(reflla.lon)*(-satpos.r*cos(satpos.e)*cos(satpos.a))
                -sin(reflla.lon)*(satpos.r*cos(satpos.e)*sin(satpos.a))
                +cos(reflla.lat)*cos(reflla.lon)*(satpos.r*sin(satpos.e)) + refpos.x;
    
    ecefret.y=  sin(reflla.lat)*sin(reflla.lon)*(-satpos.r*cos(satpos.e)*cos(satpos.a))
                +cos(reflla.lon)*(satpos.r*cos(satpos.e)*sin(satpos.a))
                +cos(reflla.lat)*sin(reflla.lon)*(satpos.r*sin(satpos.e)) + refpos.y;

    ecefret.z=  -cos(reflla.lat)*(-satpos.r*cos(satpos.e)*cos(satpos.a))
                +0
                +sin(reflla.lat)*(satpos.r*sin(satpos.e)) + refpos.z;

    return ecefret;
};


ear ECEFtoEAR(ecef satpos, lla reflla){
    // returns the EAR values seen by a ground station of a satellite at the given ECEF position
    ear earret;
    ecef refpos= LLAtoECEF(reflla);
    ecef relpos;
    relpos.x=-satpos.x+refpos.x;
    relpos.y=-satpos.y+refpos.y;
    relpos.z=-satpos.z+refpos.z;

    earret.e=   acos((refpos.x*relpos.x+refpos.y*relpos.y+refpos.z*relpos.z)
                /sqrt((pow2(refpos.x)+pow2(refpos.y)+pow2(refpos.z))*(pow2(relpos.x)+pow2(relpos.y)+pow2(relpos.z))))-M_PI_2;

    earret.a=   atan((-refpos.y*relpos.x+refpos.x*relpos.y/sqrt((pow2(refpos.x)+pow2(refpos.y)+pow2(refpos.z))*(pow2(relpos.x)+pow2(relpos.y)+pow2(relpos.z))))
                /(-refpos.z*refpos.x*relpos.x-refpos.z*refpos.y*relpos.y+(pow2(refpos.x)+pow2(refpos.y))*relpos.z)/sqrt((pow2(refpos.x)+pow2(refpos.y))+(pow2(refpos.x)+pow2(refpos.y)+pow2(refpos.z))*(pow2(relpos.x)+pow2(relpos.y)+pow2(relpos.z))));

    earret.r=   sqrt(pow2(relpos.x)+pow2(relpos.y)+pow2(relpos.z));
    return earret;
};

// thcs EARtoSEZ(ear argpos, ear argvel, eci refpos) {// need to double check this maths
//     //also probably better suited elsewhere as this is a tracking calculation
//     thcs sezposret;
//     thcs sezvelret;
//     sezposret.s=refpos.i- argpos.r*cos(argpos.e)*cos(argpos.a);
//     sezposret.e=refpos.j- argpos.r*cos(argpos.e)*sin(argpos.a);
//     sezposret.z=refpos.k- argpos.r*sin(argpos.e);

//     sezvelret.s=-argvel.r*cos(argpos.e)*cos(argpos.a)+argpos.r*sin(argpos.e)*cos(argpos.a)*argvel.e+argpos.r*cos(argpos.e)*sin(argpos.a)*argvel.a;
//     sezvelret.e=argvel.r*cos(argpos.e)*sin(argpos.a)-argpos.r*sin(argpos.e)*sin(argpos.a)*argvel.e+argpos.r*cos(argpos.e)*cos(argpos.a)*argvel.a;
//     sezvelret.z=argvel.r*sin(argpos.e)*argpos.r*cos(argpos.e)*argvel.e;
//     return sezposret;
//     return sezvelret;
// };

ecef ENUtoECEF(enu satpos, lla reflla){
    // returns the ECEF position of a satellite tracked by a ground station in East North Up
    ecef ecefret;
    ecef refpos = LLAtoECEF(reflla);

    ecefret.x=  -sin(reflla.lon)*satpos.e
                -sin(reflla.lat)*cos(reflla.lon)*satpos.n
                +cos(reflla.lat)*cos(reflla.lon)*satpos.u + refpos.x;

    ecefret.y=  cos(reflla.lon)*satpos.e
                -sin(reflla.lat)*sin(reflla.lon)*satpos.n
                +cos(reflla.lat)*sin(reflla.lon)*satpos.u + refpos.y;

    ecefret.z=  0
                +cos(reflla.lat)*satpos.n
                +sin(reflla.lat)*satpos.u + refpos.z;

    return ecefret;
};

enu ECEFtoENU(ecef satpos, lla reflla){
    // returns the ENU values seen by a ground station of a satellite at the given ECEF position
    enu enuret;
    ecef refpos = LLAtoECEF(reflla);

    enuret.e=   -sin(reflla.lon)*(satpos.x-refpos.x)
                +cos(reflla.lon)*(satpos.y-refpos.y)
                +0;

    enuret.n=   -sin(reflla.lat)*cos(reflla.lon)*(satpos.x-refpos.x)
                -sin(reflla.lat)*sin(reflla.lon)*(satpos.y-refpos.y)
                +cos(reflla.lat)*(satpos.z-refpos.z);

    enuret.u=   cos(reflla.lat)*cos(reflla.lon)+(satpos.x-refpos.x)
                +cos(reflla.lat)*sin(reflla.lon)*(satpos.y-refpos.y)
                +sin(reflla.lat)*(satpos.z-refpos.z);

    return enuret;
};

ecef SEZtoECEF(thcs satpos, lla reflla){
    // returns the ECEF position of a satellite tracked by a ground station in the Topocentric Horizon Co-ordinate System
    ecef ecefret;
    ecef refpos = LLAtoECEF(reflla);

    ecefret.x=  cos(reflla.lon)*satpos.s
                +sin(reflla.lat)*cos(reflla.lon)*satpos.e
                + refpos.z;

    ecefret.y=  -sin(reflla.lat)*sin(reflla.lon)*satpos.s
                +sin(reflla.lat)*cos(reflla.lon)*satpos.e
                +cos(reflla.lat)*satpos.z + refpos.y;

    ecefret.z=  cos(reflla.lat)*sin(reflla.lon)*satpos.s
                -cos(reflla.lat)*cos(reflla.lon)*satpos.e
                +sin(reflla.lat)*satpos.z + refpos.z;

    return ecefret;
};

thcs ECEFtoSEZ(ecef satpos, lla reflla){
    // returns the SEZ values seen by a ground station of a satellite at the given ECEF position
    thcs sezret;
    ecef refpos = LLAtoECEF(reflla);

    sezret.s=   sin(reflla.lat)*cos(reflla.lon)*(satpos.x-refpos.x)
                +sin(reflla.lat)*sin(reflla.lon)*(satpos.y-refpos.y)
                -cos(reflla.lat)*(satpos.z-refpos.z);

    sezret.e=   -sin(reflla.lon)*(satpos.x-refpos.x)
                +cos(reflla.lon)*(satpos.y-refpos.y)
                +0;

    sezret.z=   cos(reflla.lat)*cos(reflla.lon)+(satpos.x-refpos.x)
                +cos(reflla.lat)*sin(reflla.lon)*(satpos.y-refpos.y)
                +sin(reflla.lat)*(satpos.z-refpos.z);

    return sezret;
};

eci PCStoECI(orbparam arg, pcs pcsarg){
    eci eciret;

    eciret.i=   (cos(arg.asc)*cos(arg.aop)-sin(arg.asc)*sin(arg.aop)*cos(arg.inc))*pcsarg.p +
                (-cos(arg.asc)*sin(arg.aop)-cos(arg.asc)*sin(arg.aop)*cos(arg.inc))*pcsarg.q +
                (sin(arg.asc)*sin(arg.inc))*pcsarg.w;

    eciret.j=   (sin(arg.asc)*cos(arg.aop)-cos(arg.asc)*sin(arg.aop)*cos(arg.inc))*pcsarg.p +
                (-sin(arg.asc)*sin(arg.aop)-cos(arg.asc)*cos(arg.aop)*cos(arg.inc))*pcsarg.q +
                (-cos(arg.asc)*sin(arg.inc))*pcsarg.w;

    eciret.k=   sin(arg.aop)*sin(arg.inc)*pcsarg.p +
                cos(arg.aop)*sin(arg.inc)*pcsarg.q +
                cos(arg.inc)*pcsarg.w;

    eciret.vi=   (cos(arg.asc)*cos(arg.aop)-sin(arg.asc)*sin(arg.aop)*cos(arg.inc))*pcsarg.vp +
                (-cos(arg.asc)*sin(arg.aop)-cos(arg.asc)*sin(arg.aop)*cos(arg.inc))*pcsarg.vq +
                (sin(arg.asc)*sin(arg.inc))*pcsarg.vw;

    eciret.vj=   (sin(arg.asc)*cos(arg.aop)-cos(arg.asc)*sin(arg.aop)*cos(arg.inc))*pcsarg.vp +
                (-sin(arg.asc)*sin(arg.aop)-cos(arg.asc)*cos(arg.aop)*cos(arg.inc))*pcsarg.vq +
                (-cos(arg.asc)*sin(arg.inc))*pcsarg.vw;

    eciret.vk=   sin(arg.aop)*sin(arg.inc)*pcsarg.vp +
                cos(arg.aop)*sin(arg.inc)*pcsarg.vq +
                cos(arg.inc)*pcsarg.vw;

    return eciret;
};

pcs KOEtoPCS(orbparam arg){
    // gives the ECI position and velocity of a satellite given its Kepler Orbital Elements

    pcs pcsret;
    double orbitradius;
    orbitradius = (arg.sma*(1-pow2(arg.ecc)))/(1+arg.ecc*cos(arg.truanom));
    pcsret.p = orbitradius * cos(arg.truanom); 
    pcsret.q = orbitradius * sin(arg.truanom);
    pcsret.w = 0;

    double p = arg.sma*(1-pow2(arg.ecc));
    double thetadot = sqrt(planet.sgp*p)/pow2(orbitradius);
    double rdot = sqrt(planet.sgp/p)*arg.ecc*sin(arg.truanom);

    pcsret.vp = rdot*cos(arg.truanom)-orbitradius*sin(arg.truanom)*thetadot;
    pcsret.vq = rdot*sin(arg.truanom)-orbitradius*cos(arg.truanom)*thetadot;
    pcsret.vw = 0;
    return pcsret;
};

orbparam ECItoKOE(eci posvel){
    orbparam koeret;
    eci argangmnt;
    eci nodalvect;
    pcs eccvect;
    double orbitalradius = sqrt(pow2(posvel.i)+pow2(posvel.j)+pow2(posvel.k));
    double orbitalvelocity = sqrt(pow2(posvel.vi)+pow2(posvel.vj)+pow2(posvel.vk));
    argangmnt.i = (posvel.j*posvel.vk-posvel.vj*posvel.k);
    argangmnt.j = (posvel.k*posvel.vi-posvel.i*posvel.vk);
    argangmnt.k = (posvel.i*posvel.vj-posvel.j*posvel.vi);
    double argangmntnorm=sqrt(pow2(argangmnt.i)+pow2(argangmnt.j)+pow2(argangmnt.k));
    koeret.inc = acos(argangmnt.k/argangmntnorm);
    nodalvect.i=-argangmnt.j;
    nodalvect.j=argangmnt.i;
    double nodalvectnorm =  atan2(nodalvect.i,nodalvect.j);
    koeret.asc = acos(nodalvect.i/nodalvectnorm);
    //if Ny > 0 then 0<asc<180
    //if Ny < 0 then 180<asc<360
    if (nodalvect.j<0){
        koeret.asc=2*M_PI-koeret.asc;
    };
    eccvect.p = (1/planet.sgp)*pow2(orbitalvelocity);
    eccvect.q = (-1/planet.sgp)*(posvel.i*posvel.vi+posvel.j*posvel.j+posvel.k*posvel.vk);
    eccvect.w = 0;
    eccvect.vp = 0; eccvect.vq = 0; eccvect.vw;
    eci ecceci = PCStoECI(koeret, eccvect);// this may cause an error later
    koeret.ecc = sqrt(pow2(eccvect.p)+pow2(eccvect.q));
    koeret.sma=(pow2(argangmntnorm)/planet.sgp)/(1-pow2(koeret.ecc));
    koeret.aop=acos((nodalvect.i*eccvect.p+nodalvect.j*eccvect.q)/(nodalvectnorm*koeret.ecc));
    //if ecc_z > 0 then 0<aop<180
    //if ecc_z < 0 then 180<aop<360
    if(ecceci.k<0){
        koeret.aop=2*M_PI-koeret.aop;
    };
    koeret.truanom=acos((koeret.sma*(1-pow2(koeret.ecc))-orbitalradius)/(koeret.ecc*orbitalradius));
    if(posvel.i*posvel.vi+posvel.j*posvel.vj*posvel.k*posvel.vk<0){//dot product of position and velocity
        koeret.truanom=2*M_PI-koeret.truanom;
    };
    koeret.eccanom=atan2(sqrt(1-koeret.ecc)*sin(koeret.truanom/2),sqrt(1+koeret.ecc)*cos(koeret.truanom/2));
    koeret.meananom=koeret.eccanom-koeret.ecc*sin(koeret.eccanom);
    return koeret;
};

ecef ECItoECEF(eci arg,double siderealangle){
    ecef ecefret;
    ecefret.x = cos(siderealangle)*arg.i-sin(siderealangle)*arg.j;
    ecefret.y = sin(siderealangle)*arg.i+cos(siderealangle)*arg.j;
    ecefret.z = arg.k; //ignoring precession and nutation
    return ecefret;
};

eci ECEFtoECI(ecef arg, double siderealangle){
    eci eciret;

    eciret.i= cos(siderealangle)*arg.x+sin(siderealangle)*arg.y;
    eciret.j=-sin(siderealangle)*arg.x+cos(siderealangle)*arg.y;
    eciret.k=arg.z; //ignoring precession and nutation
    return eciret;
};

//the following two are not exactly true axis transforms, and are only used for maneuvers and pointing command handling

eci VNBtoECI(eci posvel, vnb VNBdV){
    eci ECIdV;
    //precalculating as this needs to be run quite fast
    eci rxv, vxrxv;
    rxv.i=posvel.j*posvel.vk-posvel.k*posvel.vj;
    rxv.j=posvel.k*posvel.vi-posvel.i*posvel.vk;
    rxv.k=posvel.i*posvel.vj-posvel.j*posvel.vi;

    vxrxv.i=posvel.vj*rxv.k-posvel.vk*rxv.j;
    vxrxv.j=posvel.vk*rxv.i-posvel.vi*rxv.k;
    vxrxv.k=posvel.vi*rxv.j-posvel.vj*rxv.i;

    double absv, absrxv, absvxrxv;
    absv= sqrt(pow2(posvel.vi)+pow2(posvel.vj)+pow2(posvel.vk));
    absrxv=sqrt(pow2(rxv.i)+pow2(rxv.j)+pow2(rxv.k));
    absvxrxv=sqrt(pow2(vxrxv.i)+pow2(vxrxv.j)+pow2(vxrxv.k));



    ECIdV.vi       =(posvel.vi/absv)*VNBdV.v
                    +(rxv.i/absrxv)*VNBdV.n
                    +(vxrxv.i/absvxrxv)*VNBdV.b;

    ECIdV.vi       =(posvel.vj/absv)*VNBdV.v
                    +(rxv.j/absrxv)*VNBdV.n
                    +(vxrxv.j/absvxrxv)*VNBdV.b;

    ECIdV.vi       =(posvel.vk/absv)*VNBdV.v
                    +(rxv.k/absrxv)*VNBdV.n
                    +(vxrxv.k/absvxrxv)*VNBdV.b;
    return ECIdV;
};

vnb ECItoVNB(eci posvel, eci ECIdV){
    vnb VNBdV;
    //precalculating as this needs to be run quite fast
    eci rxv, vxrxv;
    rxv.i=posvel.j*posvel.vk-posvel.k*posvel.vj;
    rxv.j=posvel.k*posvel.vi-posvel.i*posvel.vk;
    rxv.k=posvel.i*posvel.vj-posvel.j*posvel.vi;

    vxrxv.i=posvel.vj*rxv.k-posvel.vk*rxv.j;
    vxrxv.j=posvel.vk*rxv.i-posvel.vi*rxv.k;
    vxrxv.k=posvel.vi*rxv.j-posvel.vj*rxv.i;

    double absv, absrxv, absvxrxv;
    absv= sqrt(pow2(posvel.vi)+pow2(posvel.vj)+pow2(posvel.vk));
    absrxv=sqrt(pow2(rxv.i)+pow2(rxv.j)+pow2(rxv.k));
    absvxrxv=sqrt(pow2(vxrxv.i)+pow2(vxrxv.j)+pow2(vxrxv.k));

    VNBdV.v       =(posvel.vi/absv)*ECIdV.vi+(posvel.vj/absv)*ECIdV.vj+(posvel.vk/absv)*ECIdV.vk;

    VNBdV.n       =(rxv.i/absrxv)*ECIdV.vi+(rxv.j/absrxv)*ECIdV.vj+(rxv.k/absrxv)*ECIdV.vk;

    VNBdV.b       =(vxrxv.i/absvxrxv)*ECIdV.vi+(vxrxv.j/absvxrxv)*ECIdV.vj+(vxrxv.k/absvxrxv)*ECIdV.vk;
    return VNBdV;
};

};

    