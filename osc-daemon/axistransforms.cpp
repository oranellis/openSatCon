#include <iostream>
#include <vector>
#include <math.h>

#include "osctypes.hpp"

namespace osc{
double meantotrue(orbparam KOE){
        if (orbparam.ecc < 0.2){
            orbparam.truanon = orbparam.meananom+2*orbparam.ecc*sin(orbparam.meananom)+1.25*orbparam.ecc^2*sin(2*orbparam.meananom)-orbparam.ecc^3*(0.25*sin(orbparam.meananom)-(13/12)*sin(3*orbparam.meananom))
        } else{
            // newton-raphson's method
        };
};


double greenwichsiderealangle(utc){
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
    double e=2.0*flattening-flattening^2;
    double rho; //not sure what rho is yet
    double R;


    lat=atan(arg.z/(sqrt(arg.x^2+arg.y^2)));
    llaret.lon=atan(arg.y/arg.x)-siderealangle;
    double lat_=lat+1.0e6;
    R=rho*cos(lat);
    while (abs(lat-lat_) > 1.0e-6){
        C=1/sqrt(1-e*sin(lat)^2);
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
    orbitradius = (arg.sma*(1-arg.ecc^2))/(1+arg.ecc*cos(arg.truanom));
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
    double mu = 10; //find value
    double orbitalradius = sqrt(argpos.x^2+argpos.y^2+argpos.z^2);
    double orbitalvelocity = sqrt(argvel.x^2+argvel.y^2+argvel.z^2);
    //angular momentum vector is wrong
    argangmnt.x = (argpos.y*argvel.z-argvel.y*argpos.z);
    argangmnt.x = (argpos.z*argvel.x-argpos.x*argvel.z);
    argangmnt.x = (argpos.x*argvel.y-argpos.y*argvel.x);
    //h_norm=norm(angularmomentumvector)
    //i = acos(hz/h);
    //N_hat=cross(z_hat, h_norm)
    //asc = acos(Nx/N)
    //if Ny > 0 then 0<asc<180
    //if Ny < 0 then 180<asc<360
    //ecc = (1/mu)*((v^2-(mu/r))r_hat-dot(r,v)v_hat)
    //ecc_norm =norm(ecc)
    //sma=(h^2/mu)/(1-ecc^2)
    //aop=acos(dot(N, ecc)/N*ecc)
    //if ecc_z > 0 then 0<aop<180
    //if ecc_z < 0 then 180<aop<360
    //truanom=acos(dot(ecc, r)/e*r)
    //if dot(r, v) > 0 then 0<trueanom<180
    //if dot(r, v) < 0 then 180<trueanom<360
    //return koeret;
};
};

    