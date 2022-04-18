#ifndef OSCTYPES_H
#define OSCTYPES_H

#include <vector>
#include <array>
#include <math.h>

namespace osc {

    // Common data types
    struct vec3 {
        // Class variable initialisers
        std::array<double, 3> data = {0, 0, 0};

        // Initialisers
        vec3(std::array<double, 3> initArray):data(initArray) {}

        vec3(double i, double j, double k):data({i, j, k}) {}

        vec3() {}

        // Operator overloads
        vec3 operator+ (vec3 rhs) { return vec3(data[0] + rhs[0], data[1] + rhs[1], data[2] + rhs[2]); }

        vec3 operator- (vec3 rhs) { return vec3(data[0] - rhs[0], data[1] - rhs[1], data[2] - rhs[2]); }

        vec3 operator* (double rhs) { return vec3(data[0] * rhs, data[1] * rhs, data[2] * rhs); }

        vec3 operator/ (double rhs) { return vec3(data[0] / rhs, data[1] / rhs, data[2] / rhs); }

        double operator [] (int i) { return data[i]; }

        // Implicit type conversion
        operator std::array<double, 3> () { return data; };

        // Member functions
        vec3 dot(vec3 arg) {
            /*
            Performs the vector dot product with the argument
            */
            return vec3(data[0] * arg[0], data[1] * arg[1], data[2] * arg[2]);
        }

        vec3 cross(vec3 arg) {
            /*
            Performs the vector cross product with the argument
            */
            return vec3(data[1]*arg[2] - data[2]*arg[1], data[2]*arg[0] - data[0]*arg[2], data[0]*arg[1] - data[1]*arg[0]);
        }
    };



    struct position {
        // Class variable initialisers
        double x = 0;
        double y = 0;
        double z = 0;

        // Initialisers
        position(double initX, double initY, double initZ):x(initX), y(initY), z(initZ) {}

        position(std::array<double, 3> initPos):x(initPos[0]), y(initPos[1]), z(initPos[2]) {}
        
        position() {}

        // Operator overloads
        position operator+(vec3 rhs) {
            /*
            Addition overload for addition of two positions
            */
            return position(x + rhs[0], y + rhs[1], z + rhs[2]);
            vec3 thiun, aalsdkfj;
        }

        position operator-(vec3 rhs) {
            /*
            Subtraction overload for subtraction of two positions
            */
            return position(x - rhs[0], y - rhs[1], z - rhs[2]);
        }

        // Implicit type conversions
        operator std::array<double, 3> () { return {x, y, z}; }

        operator vec3 () { return {x, y, z}; }

        // Member functions
        position multiply(const double rhs) {
            /*
            Multiplication by a constant
            */
            return position(x * rhs, y * rhs, z * rhs);
        }

        position divide(const double rhs) {
            /*
            Division by a constant
            */
            return position(x / rhs, y / rhs, z / rhs);
        }

        position dot(vec3 argPos) {
            /*
            Performs the vector dot product with the argument
            */
            return position(x * argPos[0], y * argPos[1], z * argPos[2]);
        }
    };



    struct momentofinertia {
        // Class variable initialisers
        double Ixx;
        double Iyy;
        double Izz;

        // Member functions
        void addMass(momentofinertia I, double m, position r) {
            /*
            Adds moment of inertia of a component with added with mass m, distance from the overall CG r and moment of inertia I
            */
            Ixx += I.Ixx + m * r.y * r.y * r.z * r.z;
            Iyy += I.Iyy + m * r.z * r.z * r.x * r.x;
            Izz += I.Izz + m * r.x * r.x * r.x * r.x;
        }
    };



    struct quaternion {
        // Class variable initialisers
        double qw = 0;
        double qx = 0;
        double qy = 0;
        double qz = 0;

        // Initialisers
        quaternion() {}

        quaternion(std::array<double, 3> vec1, std::array<double, 3> vec2) {
            /*
            Generates a quaternion rotation between two unit vectors in the same reference frame
            */
            qw =    sqrt(pow(pow(vec1[0],2)+pow(vec1[1],2)+pow(vec1[2],2),2)*(pow(pow(vec2[0],2)+pow(vec2[1],2)+pow(vec2[2],2),2)));
                            +vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
            qx =    vec1[1]*vec2[2]-vec1[2]*vec2[1];
            qy =    vec1[2]*vec2[0]-vec1[0]*vec2[2];
            qz =    vec1[0]*vec2[1]-vec1[1]*vec2[0];
        }

        // Member functions
        std::array<std::array<double, 3>, 3> toMat() {
            /*
            Information for conversion available at:
            https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
            */
            std::array<std::array<double, 3>, 3> result;
            result[0] = {{1-2*qy*qy-2*qz*qz, 2*qx*qy-2*qz*qw, 2*qx*qz+2*qy*qw}};
            result[1] = {{2*qx*qy+2*qz*qw, 1-2*qx*qx-2*qz*qz, 2*qy*qz-2*qx*qw}};
            result[2] = {{2*qx*qz-2*qy*qw, 2*qy*qz+2*qx*qw, 1-2*qx*qx-2*qy*qy}};
            return result;
        }
    };
    
    
    
    struct forceTorqueModel {
        /*
        Represents the vector of forces and moments of an actuator at maximum actuation
        */
        
        // Member vars
        std::array<double, 6> ftVec {0, 0, 0, 0, 0, 0};
        
        // Initialisers
        forceTorqueModel(double Fx, double Fy, double Fz, double Txx, double Tyy, double Tzz):ftVec({Fx, Fy, Fz, Txx, Tyy, Tzz}) {
            /* 
            Generic initialiser for the vector, directly assigning each component
            */
        }

        forceTorqueModel(vec3 maxThrust, position thrusterPos) {
            /*
            Initialiser for thrusters with a thrust magnitude and direction of action
            */
            ftVec[0] = maxThrust[0]; // Resultant force on the object cg is off centre force
            ftVec[1] = maxThrust[1];
            ftVec[2] = maxThrust[2];

            vec3 torque = ((vec3)thrusterPos).cross(maxThrust);

            ftVec[3] = torque[0]; // Resultant force on the object cg is off centre force
            ftVec[4] = torque[1]; // Resultant force on the object cg is off centre force
            ftVec[5] = torque[2]; // Resultant force on the object cg is off centre force
        }

        forceTorqueModel(std::array<double, 6> initFTV):ftVec(initFTV) {}

        // Implicit type converters
        operator std::array<double, 6> () const { 
            /*
            Adds support for implicit conversion from forceTorqueModel to std::array<double, 6>
            */
            return ftVec; 
        }

        // Member functions
        forceTorqueModel normalise() {
            double mag = 1/sqrt( pow2(ftVec[0]) + pow2(ftVec[0]) + pow2(ftVec[0]) + pow2(ftVec[0]) + pow2(ftVec[0]) + pow2(ftVec[0]) );
            return std::array<double, 6> { ftVec[0]*mag, ftVec[1]*mag, ftVec[2]*mag, ftVec[3]*mag, ftVec[4]*mag, ftVec[5]*mag };
        }
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



    // Inline helper functions
    inline double pow2(double arg) { return arg*arg; }
}

#endif // OSCTYPES_H