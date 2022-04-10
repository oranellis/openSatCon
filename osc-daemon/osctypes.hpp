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

    struct quaternion {
        double qw;
        double qx;
        double qy;
        double qz;

        std::array<std::array<double, 3>, 3> toMat() {}
    };

    struct orbparam {
        long int sma; // semi major axis (m)
        double ecc; // eccentricity
        double inc; // inclination (rad)
        double aop; // argument of periapsis (rad)
        double asc; // longitude of the ascending node (rad)
        double meananom; // mean anomaly (rad)
        double eccanom; // eccentric anomaly (rad)
        double truanon; // true anomaly (rad)
        double trulon; // true longitude (rad)
    };

    struct pcspos {
        // Orbital position of satellite in Perifocal Co-Ordinate System - satellite centred
        double p; // points towards periapsis of orbit
        double q; // right hand angle along orbital plane
        double w; // normal to orbital plane
    };

    struct ecipos {
        // Orbital position of satellite in Earth Centred Inertial Co-Ordinate System
        double i; // vector from Earth to sun on J2000, 2000-01-01 at 12:00 TT
        double j; // orthogonal towards Equator
        double k; // passes through Celestial North Pole
    };

    struct ecefpos {
        // Orbital position of satellite in Earth Centred Earth Fixed Co-Ordinate System
        double x; // vector passing through Greenwich Meridian
        double y; // orthogonal towards Equator
        double z; // passes through Celestial North Pole
    };

    struct nedpos {
        // Orbital position of satellite in North East Down Co-Ordinate System - satellite centred
        double n; // points North
        double e; // points East
        double d; // ppoints down
    };

    struct thcspos {
        // Orbital position of satellite in Topocentric Horizon Co-Ordinate System - centred on ground point
        double s; // points South
        double e; // points East
        double z; // points up
    };

    struct earpos {
        // Orbital position of satellite in Elevation Azimuth Range Co-Ordinate System - centred on ground point
        double e; // elevation (rad)
        double a; // azimuth (rad)
        double r; // points up
    };

    struct orbrot {
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