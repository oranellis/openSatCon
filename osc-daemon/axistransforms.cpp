#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>

#include "osctypes.hpp"
#include "planet.cpp"

namespace osc{
double meantotrue(orbparam KOE){
        if (KOE.ecc < 0.2){
            KOE.truanom = KOE.meananom+2*KOE.ecc*sin(KOE.meananom)+1.25*pow(KOE.ecc,2)*sin(2*KOE.meananom)-pow(KOE.ecc,3)*(0.25*sin(KOE.meananom)-(13/12)*sin(3*KOE.meananom));
        } else{
            // newton-raphson's method
        };
};


double greenwichsiderealangle(){
    time_t timer;
    //for Oran
    struct tm y2k = {0};
    double seconds;

    y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
    y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;

    time(&timer);

    seconds = difftime(timer,mktime(&y2k));


};



ecef LLAtoECEF(lla arg) {
    ecef ecefret;
    // converts the angular position of the satellite to the ECEF co-ordinate system
    // this gives the correct altitude of a non-spherical earth, note that altitude here is above the ground
    double normaldistance=planet.semimajoraxis/sqrt(1-(planet.eccentricity*sin(arg.lat)));
    ecefret.x=(normaldistance+arg.alt)*cos(arg.lon)*cos(arg.lat);
    ecefret.y=(normaldistance+arg.alt)*sin(arg.lon)*cos(arg.lat);
    ecefret.z=(normaldistance*(1-pow(planet.eccentricity,2))+arg.alt)*sin(arg.lat);
    return ecefret;
};

lla ECEFtoLLA(ecef arg){
    //returns ground position of satellite from ECEF co-ords
    lla llaret;
    double secondeccentricity=planet.eccentricity/sqrt(1-pow(planet.eccentricity,2));
    double p = sqrt(pow(arg.x,2)+pow(arg.y,2));
    double theta = atan2(arg.z*planet.semimajoraxis,p*planet.semiminoraxis);
    llaret.lon = atan2(arg.y,arg.x);
    llaret.lat = atan2(arg.z+pow(secondeccentricity,2)*planet.semiminoraxis*pow(sin(theta),3),p-pow(planet.eccentricity,2)*planet.semimajoraxis*pow(sin(theta),3));
    double normaldistance = planet.semimajoraxis/sqrt(1-pow(planet.eccentricity,2)*pow(sin(llaret.lat),2));
    llaret.alt = (p/cos(llaret.lat))-normaldistance;
    return llaret;
};

ned ECEFtoNED(ecef satpos, ecef refpos){ //probably better to use LLA for reference position
    ned nedret;
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

ecef NEDtoECEF(ned satpos, lla refpos){
    ecef ecefret;

    ecefret.x=  -cos(refpos.lon)*sin(refpos.lat)*satpos.n
                -sin(refpos.lon)*satpos.e
                -cos(refpos.lon)*cos(refpos.lat)*satpos.d;

    ecefret.y=  -sin(refpos.lon)*sin(refpos.lat)*satpos.n
                +cos(refpos.lon)*satpos.e
                -sin(refpos.lon)*cos(refpos.lat)*satpos.d;

    ecefret.z=  cos(refpos.lat)*satpos.n
                +0
                -sin(refpos.lat)*satpos.d;
    return ecefret;   
};

ecef EARtoECEF(ear satpos, lla refpos){//this method goes through SEZ co-ordinates, but seems to be the best method
    ecef ecefret;
    ecef ecefrefpos = LLAtoECEF(refpos);
    ecefret.x=  sin(refpos.lat)*cos(refpos.lon)*(-satpos.r*cos(satpos.e)*cos(satpos.a))
                -sin(refpos.lon)*(satpos.r*cos(satpos.e)*sin(satpos.a))
                +cos(refpos.lat)*cos(refpos.lon)*(satpos.r*sin(satpos.e)) + ecefrefpos.x;
    
    ecefret.y=  sin(refpos.lat)*sin(refpos.lon)*(-satpos.r*cos(satpos.e)*cos(satpos.a))
                +cos(refpos.lon)*(satpos.r*cos(satpos.e)*sin(satpos.a))
                +cos(refpos.lat)*sin(refpos.lon)*(satpos.r*sin(satpos.e)) + ecefrefpos.y;

    ecefret.z=  -cos(refpos.lat)*(-satpos.r*cos(satpos.e)*cos(satpos.a))
                +0
                +sin(refpos.lat)*(satpos.r*sin(satpos.e)) + ecefrefpos.z;

    return ecefret;
};


ear ECEFtoEAR(ecef satpos, ecef refpos){ //probably better to use LLA for reference position
    ear earret;
    ecef relpos;
    relpos.x=-satpos.x+refpos.x;
    relpos.y=-satpos.y+refpos.y;
    relpos.z=-satpos.z+refpos.z;

    earret.e=   acos((refpos.x*relpos.x+refpos.y*relpos.y+refpos.z*relpos.z)
                /sqrt((pow(refpos.x,2)+pow(refpos.y,2)+pow(refpos.z,2))*(pow(relpos.x,2)+pow(relpos.y,2)+pow(relpos.z,2))))-M_PI_2;

    earret.a=   atan((-refpos.y*relpos.x+refpos.x*relpos.y/sqrt((pow(refpos.x,2)+pow(refpos.y,2)+pow(refpos.z,2))*(pow(relpos.x,2)+pow(relpos.y,2)+pow(relpos.z,2))))
                /(-refpos.z*refpos.x*relpos.x-refpos.z*refpos.y*relpos.y+(pow(refpos.x,2)+pow(refpos.y,2))*relpos.z)/sqrt((pow(refpos.x,2)+pow(refpos.y,2))+(pow(refpos.x,2)+pow(refpos.y,2)+pow(refpos.z,2))*(pow(relpos.x,2)+pow(relpos.y,2)+pow(relpos.z,2))));

    earret.r=   sqrt(pow(relpos.x,2)+pow(relpos.y,2)+pow(relpos.z,2));
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

ecef ENUtoECEF(enu satpos, lla refpos){
    ecef ecefret;
    ecef ecefrefpos = LLAtoECEF(refpos);

    ecefret.x=  -sin(refpos.lon)*satpos.e
                -sin(refpos.lat)*cos(refpos.lon)*satpos.n
                +cos(refpos.lat)*cos(refpos.lon)*satpos.u + ecefrefpos.x;

    ecefret.y=  cos(refpos.lon)*satpos.e
                -sin(refpos.lat)*sin(refpos.lon)*satpos.n
                +cos(refpos.lat)*sin(refpos.lon)*satpos.u + ecefrefpos.y;

    ecefret.z=  0
                +cos(refpos.lat)*satpos.n
                +sin(refpos.lat)*satpos.u + ecefrefpos.z;

    return ecefret;
};

enu ECEFtoENU(ecef satpos, lla refpos){
    enu enuret;
    ecef ecefrefpos = LLAtoECEF(refpos);

    enuret.e=   -sin(refpos.lon)*(satpos.x-ecefrefpos.x)
                +cos(refpos.lon)*(satpos.y-ecefrefpos.y)
                +0;

    enuret.n=   -sin(refpos.lat)*cos(refpos.lon)*(satpos.x-ecefrefpos.x)
                -sin(refpos.lat)*sin(refpos.lon)*(satpos.y-ecefrefpos.y)
                +cos(refpos.lat)*(satpos.z-ecefrefpos.z);

    enuret.u=   cos(refpos.lat)*cos(refpos.lon)+(satpos.x-ecefrefpos.x)
                +cos(refpos.lat)*sin(refpos.lon)*(satpos.y-ecefrefpos.y)
                +sin(refpos.lat)*(satpos.z-ecefrefpos.z);

    return enuret;
};

ecef SEZtoECEF(thcs satpos, lla refpos){
    ecef ecefret;
    ecef ecefrefpos = LLAtoECEF(refpos);

    ecefret.x=  cos(refpos.lon)*satpos.s
                +sin(refpos.lat)*cos(refpos.lon)*satpos.e
                + ecefrefpos.z;

    ecefret.y=  -sin(refpos.lat)*sin(refpos.lon)*satpos.s
                +sin(refpos.lat)*cos(refpos.lon)*satpos.e
                +cos(refpos.lat)*satpos.z + ecefrefpos.y;

    ecefret.z=  cos(refpos.lat)*sin(refpos.lon)*satpos.s
                -cos(refpos.lat)*cos(refpos.lon)*satpos.e
                +sin(refpos.lat)*satpos.z + ecefrefpos.z;

    return ecefret;
};

thcs ECEFtoSEZ(ecef satpos, lla refpos){
    thcs sezret;
    ecef ecefrefpos = LLAtoECEF(refpos);

    sezret.s=   sin(refpos.lat)*cos(refpos.lon)*(satpos.x-ecefrefpos.x)
                +sin(refpos.lat)*sin(refpos.lon)*(satpos.y-ecefrefpos.y)
                -cos(refpos.lat)*(satpos.z-ecefrefpos.z);

    sezret.e=   -sin(refpos.lon)*(satpos.x-ecefrefpos.x)
                +cos(refpos.lon)*(satpos.y-ecefrefpos.y)
                +0;

    sezret.z=   cos(refpos.lat)*cos(refpos.lon)+(satpos.x-ecefrefpos.x)
                +cos(refpos.lat)*sin(refpos.lon)*(satpos.y-ecefrefpos.y)
                +sin(refpos.lat)*(satpos.z-ecefrefpos.z);

    return sezret;
};

eci KOEtoECI(orbparam arg){
    eci eciret;
    pcs pcs_;
    double orbitradius;
    orbitradius = (arg.sma*(1-pow(arg.ecc,2)))/(1+arg.ecc*cos(arg.truanom));
    pcs_.p = orbitradius * cos(arg.truanom); 
    pcs_.q = orbitradius * sin(arg.truanom);
    pcs_.w = 0;

    eciret.i=   (cos(arg.asc)*cos(arg.aop)-sin(arg.asc)*sin(arg.aop)*cos(arg.inc))*pcs_.p +
                (-cos(arg.asc)*sin(arg.aop)-cos(arg.asc)*sin(arg.aop)*cos(arg.inc))*pcs_.q +
                (sin(arg.asc)*sin(arg.inc))*pcs_.w;

    eciret.j=   (sin(arg.asc)*cos(arg.aop)-cos(arg.asc)*sin(arg.aop)*cos(arg.inc))*pcs_.p +
                (-sin(arg.asc)*sin(arg.aop)-cos(arg.asc)*cos(arg.aop)*cos(arg.inc))*pcs_.q +
                (-cos(arg.asc)*sin(arg.inc))*pcs_.w;

    eciret.k=   sin(arg.aop)*sin(arg.inc)*pcs_.p +
                cos(arg.aop)*sin(arg.inc)*pcs_.q +
                cos(arg.inc)*pcs_.w;

    return eciret;
};

orbparam ECItoKOE(eci argpos, eci argvel){
    orbparam koeret;
    eci argangmnt;
    eci nodalvect;
    eccentricityvector eccvect;
    double orbitalradius = sqrt(pow(argpos.i,2)+pow(argpos.j,2)+pow(argpos.k,2));
    double orbitalvelocity = sqrt(pow(argvel.i,2)+pow(argvel.j,2)+pow(argvel.k,2));
    argangmnt.i = (argpos.j*argvel.k-argvel.j*argpos.k);
    argangmnt.j = (argpos.k*argvel.i-argpos.i*argvel.k);
    argangmnt.k = (argpos.i*argvel.j-argpos.j*argvel.i);
    double argangmntnorm=sqrt(pow(argangmnt.i,2)+pow(argangmnt.j,2)+pow(argangmnt.k,2));
    koeret.inc = acos(argangmnt.k/argangmntnorm);
    nodalvect.i=-argangmnt.j;
    nodalvect.j=argangmnt.i;
    double nodalvectnorm = sqrt(pow(nodalvect.i,2)+pow(nodalvect.j,2)); //might be atan2(nodalvect.i,nodalvect.j)
    koeret.asc = acos(nodalvect.i/nodalvectnorm);
    //if Ny > 0 then 0<asc<180
    //if Ny < 0 then 180<asc<360
    eccvect.r = (1/planet.stdgravparam)*pow(orbitalvelocity,2);
    eccvect.v = (-1/planet.stdgravparam)*(argpos.i*argvel.i+argpos.j*argpos.j+argpos.k*argvel.k);
    koeret.ecc = sqrt(pow(eccvect.r,2)+pow(eccvect.v,2));
    koeret.sma=(pow(argangmntnorm,2)/planet.stdgravparam)/(1-pow(koeret.ecc,2));
    koeret.aop=acos((nodalvect.i*eccvect.r+nodalvect.j*eccvect.v)/(nodalvectnorm*koeret.ecc));
    //if ecc_z > 0 then 0<aop<180
    //if ecc_z < 0 then 180<aop<360
    koeret.truanom=acos((koeret.sma*(1-pow(koeret.ecc,2))-orbitalradius)/(koeret.ecc*orbitalradius));
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

};

    