#include "osctypes.hpp"

namespace osc {

    // Position
    position::position(double initX, double initY, double initZ):x(initX), y(initY), z(initZ) {}

    position::position(std::array<double, 3> initPos) {
        x = initPos[0];
        y = initPos[1];
        z = initPos[2];
    }

    position position::operator+(const position& rhs) {
        /*
        Addition overload for addition of two positions
        */
        position returnPos(0, 0, 0);
        returnPos.x = x + rhs.x;
        returnPos.y = y + rhs.y;
        returnPos.z = z + rhs.z;
        return returnPos;
    }

    position position::operator-(const position& rhs) {
        /*
        Subtraction overload for subtraction of two positions
        */
        position returnPos = position();
        returnPos.x = x - rhs.x;
        returnPos.y = y - rhs.y;
        returnPos.z = z - rhs.z;
        return returnPos;
    }

    position position::multiply(const double rhs) {
        /*
        Multiplication by a constant
        */
        position returnPos = position();
        returnPos.x = x * rhs;
        returnPos.y = y * rhs;
        returnPos.z = z * rhs;
        return returnPos;
    }

    position position::divide(const double rhs) {
        /*
        Division by a constant
        */
        position returnPos = position();
        returnPos.x = x / rhs;
        returnPos.y = y / rhs;
        returnPos.z = z / rhs;
        return returnPos;
    }

    position position::dot(position argPos) {
        position returnPos = position();
        returnPos.x = x * argPos.x;
        returnPos.y = y * argPos.y;
        returnPos.z = z * argPos.z;
        return returnPos;
    }



    // momentofinertia
    momentofinertia momentofinertia::addMass(momentofinertia I, double m, position r) {
        /*
        Adds moment of inertia of a component with added with mass m, distance from the overall CG r and moment of inertia I
        */
        Ixx += I.Ixx + m * r.y * r.y * r.z * r.z;
        Iyy += I.Iyy + m * r.z * r.z * r.x * r.x;
        Izz += I.Izz + m * r.x * r.x * r.x * r.x;
    }



    // quaternion
    std::array<std::array<double, 3>, 3> quaternion::toMat() {
        /*
        Information for conversion available at:
        https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
        */
        std::array<std::array<double, 3>, 3> result;
        result[0] = {{1-2*qy*qy-2*qz*qz, 2*qx*qy-2*qz*qw, 2*qx*qz+2*qy*qw}};
        result[1] = {{2*qx*qy+2*qz*qw, 1-2*qx*qx-2*qz*qz, 2*qy*qz-2*qx*qw}};
        result[2] = {{2*qx*qz-2*qy*qw, 2*qy*qz+2*qx*qw, 1-2*qx*qx-2*qy*qy}};
    }
} // namespace osc