#include <iostream>
#include <vector>
#include <math.h>

#include "osctypes.hpp"

namespace osc{
double meantotrue(orbparam KOE){
        if (KOE.ecc < 0.2){
            KOE.truanom = KOE.meananom+2*KOE.ecc*sin(KOE.meananom)+1.25*pow(KOE.ecc,2)*sin(2*KOE.meananom)-pow(KOE.ecc,3)*(0.25*sin(KOE.meananom)-(13/12)*sin(3*KOE.meananom));
        } else{
            // newton-raphson's method
        };
};


double greenwichsiderealangle(double utc){
    //find time since J2000
    //for Oran
    double siderealangle;
    return siderealangle;
};


eci LLAtoECI(lla arg) {
    eci eciret;
    // converts the angular position of the satellite to the ECI co-ordinate system
    // note that this is from the Earth's centre and therefore geodetic latitude and geocentric longitude are not required
    eciret.i=arg.alt*cos(arg.lon);
    eciret.j=arg.alt*sin(arg.lon);
    eciret.k=arg.alt*cos(arg.lat);
    return eciret;
};

eci SEZtoECI(thcs arg) {
    //will implement later
};

lla ECEFtoLLA(ecef arg, double siderealangle){
    //returns ground position of satellite from ECEF co-ords
    lla llaret;
    double lat_=999;
    double C;
    double lat;
    double earthradius = 6378.137; 
    double flattening = 1/298.26;
    double e=(2*flattening)-(pow(flattening, 2.0));
    double rho; //not sure what rho is yet
    double R;


    lat=atan(arg.z/(sqrt(pow(arg.x,2)+pow(arg.y,2))));
    llaret.lon=atan(arg.y/arg.x)-siderealangle;
    double lat_=lat+1.0e6;
    R=rho*cos(lat);
    while (abs(lat-lat_) > 1.0e-6){
        C=1/sqrt(1-e*pow(sin(lat),2));
        lat_=atan(arg.z+earthradius*C*e*sin(lat)/R);
    }
    llaret.lat=lat_;
    llaret.alt=(R/cos(lat_)-earthradius*C);
    return llaret;
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

ecef ECItoECEF(eci arg,double siderealangle){
    ecef ecefret;
    ecefret.x = cos(siderealangle)*arg.i-sin(siderealangle)*arg.j;
    ecefret.y = sin(siderealangle)*arg.i+cos(siderealangle)*arg.j;
    ecefret.z = arg.k; //add precession and nutation
    return ecefret;
};

orbparam ECEFtoKOE(ecef argpos, ecef argvel){
    orbparam koeret;
    ecef argangmnt;
    ecef nodalvect;
    eccentricityvector eccvect;
    double mu = 10; //find value
    double orbitalradius = sqrt(pow(argpos.x,2)+pow(argpos.y,2)+pow(argpos.z,2));
    double orbitalvelocity = sqrt(pow(argvel.x,2)+pow(argvel.y,2)+pow(argvel.z,2));
    //angular momentum vector is wrong
    argangmnt.x = (argpos.y*argvel.z-argvel.y*argpos.z);
    argangmnt.y = (argpos.z*argvel.x-argpos.x*argvel.z);
    argangmnt.z = (argpos.x*argvel.y-argpos.y*argvel.x);
    double argangmntnorm=sqrt(pow(argangmnt.x,2)+pow(argangmnt.y,2)+pow(argangmnt.z,2));
    koeret.inc = acos(argangmnt.z/argangmntnorm);
    nodalvect.x=-argangmnt.y;
    nodalvect.y=argangmnt.x;
    double nodalvectnorm = sqrt(pow(nodalvect.x,2)+pow(nodalvect.y,2)); //might be atan2(nodalvect.x,nodalvect.y)
    koeret.asc = acos(nodalvect.x/nodalvectnorm);
    //if Ny > 0 then 0<asc<180
    //if Ny < 0 then 180<asc<360
    eccvect.r = (1/mu)*pow(orbitalvelocity,2);
    eccvect.v = (-1/mu)*(argpos.x*argvel.x+argpos.y*argpos.y+argpos.z*argvel.z);
    koeret.ecc = sqrt(pow(eccvect.r,2)+pow(eccvect.v,2));
    koeret.sma=(pow(argangmntnorm,2)/mu)/(1-pow(koeret.ecc,2));
    koeret.aop=acos((nodalvect.x*eccvect.r+nodalvect.y*eccvect.v)/(nodalvectnorm*koeret.ecc));
    //if ecc_z > 0 then 0<aop<180
    //if ecc_z < 0 then 180<aop<360
    double meanmotion=sqrt(mu/pow(koeret.sma,3));
    koeret.eccanom=acos(koeret.sma-orbitalradius)/(koeret.sma*orbitalradius);
    //koeret.truanom=acos(dot(ecc, r)/e*r);
        //true anomaly is a hassle
    //if dot(r, v) > 0 then 0<trueanom<180
    //if dot(r, v) < 0 then 180<trueanom<360
    return koeret;
};
};

    