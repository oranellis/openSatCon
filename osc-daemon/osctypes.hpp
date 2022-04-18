#ifndef OSCTYPES
#define OSCTYPES

#include <vector>
#include <array>

namespace osc {

    enum componentType {
    standard,
    fuelTank,
    thruster,
    rotator
    };

    struct position {
        double x;
        double y;
        double z;
        friend position operator+(const position lhs, const position& rhs) {}
        friend position operator-(const position lhs, const position& rhs) {}
        position divide(const double rhs) {}
        position multiply(const double rhs) {}
        position dot(position argPos) {}
    };

    struct momentofinertia {
        double Ix;
        double Iy;
        double Iz;
        momentofinertia addMass(momentofinertia I, double m, position r) {}
    };

    struct inertialframe {//used for attitude dynamics
        double x;
        double y;
        double z;
    };

    struct bodyframe{
        double x;
        double y;
        double z;
    };

    struct inertiamatrix {
        double  J11, J12, J13,
                J21, J22, J23, 
                J31, J32, J33;

    };

    struct quaternion {
        double qw;
        double qx;
        double qy;
        double qz;

        std::array<std::array<double, 3>, 3> toMat() {}
    };

    struct rotstates{
        double o1; //omega_i1
        double o2; //omega_i2
        double o3; //omega_i3
        double q1; //quaternion vector
        double q2;
        double q3;
    };

    struct posstates{
        double i; //eci position (m)
        double j;
        double k; 
        double vi; //eci velocity (m/s)
        double vj;
        double vk;
        double m; //spacecraft mass (kg)
    };

    struct orbparam {
        long int sma; // semi major axis (m)
        double ecc; // eccentricity
        double inc; // inclination (rad)
        double aop; // argument of periapsis (rad)
        double asc; // longitude of the ascending node (rad)
        double meananom; // mean anomaly (rad)
        double eccanom; // eccentric anomaly (rad)
        double truanom; // true anomaly (rad)
    };

    struct eccentricityvector { //remove if unnecessary
        double r; //radial component
        double v; //velocity component
    };

    struct pcs {
        // Orbital position of satellite in Perifocal Co-Ordinate System - satellite centred
        double p; // points towards periapsis of orbit
        double q; // right hand angle along orbital plane
        double w; // normal to orbital plane
        double vp, vq, vw;
    };

    struct eci {
        // Orbital position of satellite in Earth Centred Inertial Co-Ordinate System
        double i; // vector from Earth to sun on J2000, 2000-01-01 at 12:00 TT
        double j; // orthogonal towards Equator
        double k; // passes through Celestial North Pole
        double vi, vj, vk;
    };

    struct ecef {
        // Orbital position of satellite in Earth Centred Earth Fixed Co-Ordinate System
        double x; // vector passing through Greenwich Meridian
        double y; // orthogonal towards Equator
        double z; // passes through Celestial North Pole
        double vx, vy, vz;
    };

    struct ned {
        // Orbital position of satellite in North East Down Co-Ordinate System - satellite centred
        double n; // points North
        double e; // points East
        double d; // points down
        double vn, ve, vd;
    };

    struct enu {
        // Orbital position of satellite in East North Up Co-Ordinate System - Earth centred
        double e; // points East
        double n; // points North
        double u; // points Up
        double ve, vn, vu;
    };

    struct thcs {
        // Orbital position of satellite in Topocentric Horizon Co-Ordinate System - centred on ground point
        double s; // points South
        double e; // points East
        double z; // points up
        double vs, ve, vz;
    };
    
    struct lla {
        // Ground Sub-Vehicle Point on Earth's surface
        double lon; //geocentric longitude (output in deg-min-sec?)
        double lat; //geocentric latitude
        double alt; //height above WGS84 ellipsoid
    };

    struct ear {
        // Orbital position of satellite in Elevation Azimuth Range Co-Ordinate System - centred on ground point
        double e; // elevation (rad)
        double a; // azimuth (rad)
        double r; // points up
    };

        struct vnb {
        // Orbital velocity of satellite in Velocity Normal Bi-Normal axes - used for orbital maneuvers
        double v; // velocity vector - prograde
        double n; // normal vector - normal to orbital plane
        double b; // bi-normal vector - orthangonal to x and y, in the direction of the angular momentum vector
    };

    struct orbrot {//this needs fixed
        // Orbital reference rotation is using the Velocity Normal Bi-Normal system
        double q;
        double x; // velocity vector - prograde
        double y; // normal vector - normal to orbital plane
        double z; // bi-normal vector - orthangonal to x and y, in the direction of the angular momentum vector
    };

    struct orbcommand {
        double a; //altitude of apoapsis (m)
        double p; //altitude of periapsis (m)
        double e; //eccentricity
        double i; //inclination (rad)
    };

    struct powermodel {
        // Using vec indices to indicate state, 0 for low power/off, 1 for idle, 2 for in use/max load, >2 custom. Usage in W
        std::vector<double> pstates; 
    };

double pow2(double x){return x*x;}      // use this to replace pow(x,n) if n=2 and speed must =fast.lots
double pow3(double x){return x*x*x;}    // use this to replace pow(x,n) if n=3 and speed must =fast.lots
}

#endif