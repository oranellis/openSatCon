#include "osctypes.hpp"

struct position {
    double x = 0;
    double y = 0;
    double z = 0;

    friend position operator+(position lhs, const position& rhs) {
        /*
        Addition overload for addition of two positions
        */
        position returnPos = position();
        returnPos.x = lhs.x + rhs.x;
        returnPos.y = lhs.y + rhs.y;
        returnPos.z = lhs.z + rhs.z;
        return returnPos;
    }

    friend position operator-(position lhs, const position& rhs) {
        /*
        Subtraction overload for subtraction of two positions
        */
        position returnPos = position();
        returnPos.x = lhs.x - rhs.x;
        returnPos.y = lhs.y - rhs.y;
        returnPos.z = lhs.z - rhs.z;
        return returnPos;
    }

    position multiply(const double rhs) {
        /*
        Multiplication by a constant
        */
        position returnPos = position();
        returnPos.x = x * rhs;
        returnPos.y = y * rhs;
        returnPos.z = z * rhs;
        return returnPos;
    }

    position divide(const double rhs) {
        /*
        Division by a constant
        */
        position returnPos = position();
        returnPos.x = x / rhs;
        returnPos.y = y / rhs;
        returnPos.z = z / rhs;
        return returnPos;
    }

    position dot(position argPos) {
        position returnPos = position();
        returnPos.x = x * argPos.x;
        returnPos.y = y * argPos.y;
        returnPos.z = z * argPos.z;
        return returnPos;
    }
};

struct momentofinertia {
    double Ix = 0;
    double Iy = 0;
    double Iz = 0;

    momentofinertia addMass(momentofinertia I, double m, position r) {
        /*
        Adds moment of inertia of a component with added with mass m, distance from the overall CG r and moment of inertia I
        */
        Ix += I.Ix + m * r.x * r.x;
        Iy += I.Iy + m * r.y * r.y;
        Iz += I.Iz + m * r.z * r.z;
    }
};