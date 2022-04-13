#include <json/value.h>
#include <fstream>
#include <iostream>
#include <math.h>

std::ifstream craft_file("craft.json", std::ifstream::binary);
craft_file >> craft;

// access data via craft["ThrusterComponents"]["normal"] etc

//calculates moments of inertia for craft
void CalculateInertia() {
    double mass;
    double x;
    double y;
    double z;
    double Ixxtotal;
    double Iyytotal;
    double Izztotal;
    double Ixytotal;
    double Ixztotal;
    double Iyztotal;

    for(int i = 0; i < (craft["MassComponentNo"] - 1); i++)
    {
        mass = craft["MassComponents"][to_string(i)]["mass"];
        x = craft["MassComponents"][to_string(i)]["position"][0u];
        y = craft["MassComponents"][to_string(i)]["position"][1u];
        z = craft["MassComponents"][to_string(i)]["position"][2u];
        // mass calculations!
        // Ixx = m(y^2 + z^2)
        // Iyy = m(x^2 + z^2)
        // Izz = m(x^2 + y^2)
        // Ixy = -mxy
        // Ixz = -mxz
        // Iyz = -myz
        // Values are sum of moments of inertia for the point masses

        double Ixx = mass*(y^2 + z^2);
        Ixxtotal += Ixx;
        double Iyy = mass*(x^2 + z^2);
        Iyytotal += Iyy;
        double Izz = mass*(x^2 + y^2);
        Izztotal += Izz;
        double Ixy = -mass*x*y;
        Ixytotal += Ixy;
        double Ixz = -mass*x*z;
        Ixztotal += Ixz;
        double Iyz = -mass*y*z;
        Iyztotal += Iyz;
    }
}