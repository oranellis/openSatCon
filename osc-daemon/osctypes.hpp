#ifndef OSCTYPES
#define OSCTYPES

#include <vector>
#include <array>
#include <math.h>

namespace osc {

    struct position {
        double x;
        double y;
        double z;
        position(double initX, double initY, double initZ);
        position(std::array<double, 3> initPos);
        position();
        position operator+(const position& rhs);
        position operator-(const position& rhs);
        position divide(const double rhs);
        position multiply(const double rhs);
        position dot(position argPos);
        std::array<double, 3> toArray();
    };

    struct momentofinertia {
        double Ixx;
        double Iyy;
        double Izz;

        void addMass(momentofinertia I, double m, position r);
    };

    struct quaternion {
        double qw;
        double qx;
        double qy;
        double qz;

        quaternion();
        quaternion(std::array<double, 3> vec1, std::array<double, 3> vec2);
        std::array<std::array<double, 3>, 3> toMat();
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
    };

    struct eci {
        // Orbital position of satellite in Earth Centred Inertial Co-Ordinate System
        double i; // vector from Earth to sun on J2000, 2000-01-01 at 12:00 TT
        double j; // orthogonal towards Equator
        double k; // passes through Celestial North Pole
    };

    struct ecef {
        // Orbital position of satellite in Earth Centred Earth Fixed Co-Ordinate System
        double x; // vector passing through Greenwich Meridian
        double y; // orthogonal towards Equator
        double z; // passes through Celestial North Pole
    };

    struct ned {
        // Orbital position of satellite in North East Down Co-Ordinate System - satellite centred
        double n; // points North
        double e; // points East
        double d; // points down
    };

    struct enu {
        // Orbital position of satellite in East North Up Co-Ordinate System - Earth centred
        double e; // points East
        double n; // points North
        double u; // points Up
    };

    struct thcs {
        // Orbital position of satellite in Topocentric Horizon Co-Ordinate System - centred on ground point
        double s; // points South
        double e; // points East
        double z; // points up
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

    struct powermodel {
        // Using vec indices to indicate state, 0 for low power/off, 1 for idle, 2 for in use/max load, >2 custom. Usage in W
        std::vector<double> pstates; 
    };

}

#endif