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

        //vec3 operatorplus2 (vec3 term1, vec3 term2) { return vec3(data[0] + term1[0] + term2[0], data[1] + term1[1] + term2[1], data[2] + term1[2] + term2[2]); }

        vec3 operator- (vec3 rhs) { return vec3(data[0] - rhs[0], data[1] - rhs[1], data[2] - rhs[2]); }

        vec3 operator* (double rhs) { return vec3(data[0] * rhs, data[1] * rhs, data[2] * rhs); }

        vec3 operator/ (double rhs) { return vec3(data[0] / rhs, data[1] / rhs, data[2] / rhs); }

        double operator [] (int i) { return data[i]; }

        // Implicit type conversion
        operator std::array<double, 3> () { return data; };

        // Member functions
        double dot(vec3 arg) {
            /*
            Performs the vector dot product with the argument
            */
            return (data[0] * arg[0] + data[1] * arg[1] + data[2] * arg[2]);
        }

        vec3 cross(vec3 arg) {
            /*
            Performs the vector cross product with the argument
            */
            return vec3(data[1]*arg[2] - data[2]*arg[1], data[2]*arg[0] - data[0]*arg[2], data[0]*arg[1] - data[1]*arg[0]);
        }

        double mag() {
            return sqrt(dot(data));
        }

        vec3 unit() {
            return vec3(data) / vec3(data).mag();
        };
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

        quaternion(vec3 arg1, vec3 arg2) {
            vec3 cross = arg1.cross(arg2);
            qw =    sqrt(arg1.mag()*arg2.mag())+arg1.dot(arg2);
            qx =    cross[0];
            qy =    cross[1];
            qz =    cross[2];
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

        quaternion quaternionDerivative(quaternion argQuat, vec3 argRate){
            quaternion dotQuat;

            argQuat.qw=sqrt(1 - pow2(argQuat.qx) - pow2(argQuat.qy) - pow2(argQuat.qz));
            
            dotQuat.qx = 0.5 * (argQuat.qw * argRate[0]
                               -argQuat.qz * argRate[1]
                               +argQuat.qy * argRate[2]);

            dotQuat.qy = 0.5 * (argQuat.qz * argRate[0]
                               +argQuat.qw * argRate[1]
                               -argQuat.qx * argRate[2]);

            dotQuat.qz = 0.5 * (-argQuat.qy * argRate[0]
                                +argQuat.qx * argRate[1]
                                -argQuat.qw * argRate[2]);

            return dotQuat;
        };
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

    struct rotStates {
        vec3 omega; //body rates
        vec3 q; //quaternion vector part
    };

    struct posStates {
        vec3 r; //eci position (m)
        vec3 v; //eci velocity (m/s)
        double m; //spacecraft mass (kg)
    };

    struct orbParam {
        long int sma; // semi major axis (m)
        double ecc; // eccentricity
        double inc; // inclination (rad)
        double aop; // argument of periapsis (rad)
        double asc; // longitude of the ascending node (rad)
        double meanAnom; // mean anomaly (rad)
        double eccAnom; // eccentric anomaly (rad)
        double truAnom; // true anomaly (rad)
    };

    struct pcs {
        // Orbital position of satellite in Perifocal Co-Ordinate System - satellite centred
        // p - points towards periapsis of orbit
        // q - right hand angle along orbital plane
        // w - normal to orbital plane
        vec3 rPCS;
        vec3 vPCS;
    };



    struct eci {
        // Orbital position of satellite in Earth Centred Inertial Co-Ordinate System
        // i - vector from Earth to sun on J2000, 2000-01-01 at 12:00 TT
        // j - orthogonal towards Equator
        // k - passes through Celestial North Pole
        vec3 rIJK;
        vec3 vIJK;
    };



    struct ecef {
        // Orbital position of satellite in Earth Centred Earth Fixed Co-Ordinate System
        // x - vector passing through Greenwich Meridian
        // y - orthogonal towards Equator
        // z - passes through Celestial North Pole
        vec3 rXYZ;
        vec3 vXYZ;
    };



    struct ned {
        // Orbital position of satellite in North East Down Co-Ordinate System - satellite centred
        // n - points North
        // e - points East
        // d - points down
        vec3 rNED;
        vec3 vNED;
    };

    struct enu {
        // Orbital position of satellite in East North Up Co-Ordinate System - Earth centred
        // e - points East
        // n - points North
        // u - points Up
        vec3 rENU;
        vec3 vENU;
    };

    struct thcs {
        // Orbital position of satellite in Topocentric Horizon Co-Ordinate System - centred on ground point
        // s - points South
        // e - points East
        // z - points up
        vec3 rSEZ;
        vec3 vSEZ;
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
        // v - velocity vector - prograde
        // n - normal vector - normal to orbital plane
        // b - bi-normal vector - orthangonal to x and y, in the direction of the angular momentum vector
        vec3 vVNB;
    };

    struct powermodel {
        // Using vec indices to indicate state, 0 for low power/off, 1 for idle, 2 for in use/max load, >2 custom. Usage in W
        std::vector<double> pstates; 
    };
    
    // Inline helper functions
    inline double pow2(double arg) { return arg * arg; }
    inline double pow3(double arg) { return arg * arg * arg; }
}

#endif // OSCTYPES_H