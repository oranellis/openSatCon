#ifndef JSONPARSER_H
#define JSONPARSER_H

#include <iostream>
#include <list>
#include <vector>
#include <fstream>

#include "../../includes/json.hpp"
#include "../components/component.hpp"
#include "../components/fueltank.hpp"
#include "../components/actuators/thruster.hpp"
#include "../components/actuators/rotator.hpp"
#include "../osctypes.hpp"

using namespace std;
using json = nlohmann::json;

namespace osc {
    /** \mainpage openSatCon's project documentation */


    /** \struct craftconfig
    \brief craft configuration
    craftconfiguration type containing components, fueltanks, thrusters and rotators
    */
struct craftconfig {
    std::map<std::string, component> components;
    std::map<std::string, fueltank> fueltanks;
    std::map<std::string, thruster> thrusters;
    std::map<std::string, rotator> rotators;

    bool populated() {
        if (components.empty()) return false;
        else return true;
    }
};

    /** \fn parseJson(jsonPath)
    \brief parses a file at a path
    @param[in] jsonPath path to json file
    parses the json file from the path and produces a craft config object with mapped 
    components
    */
inline craftconfig parseJson(std::string jsonPath) {
    // std::cout << "start"; // debug
    // map<string, std::list<double> > Components;
    
    //json j = json::parse(jsonPath);
    std::ifstream ifs(jsonPath);
    json j = json::parse(ifs);

    craftconfig returnConfig;

    //parse and map mass components
    for(int i = 0; i < j["MassComponentNo"]; i++)
    {
        double mass = j["MassComponents"][std::to_string(i)]["mass"];
        
        momentofinertia moi;
        position pos;
        
        moi.Ixx = j["MassComponents"][std::to_string(i)]["moi"][0u];
        moi.Iyy = j["MassComponents"][std::to_string(i)]["moi"][1u];
        moi.Izz = j["MassComponents"][std::to_string(i)]["moi"][2u];
        
        pos.x = j["MassComponents"][std::to_string(i)]["position"][0u];
        pos.y = j["MassComponents"][std::to_string(i)]["position"][1u];
        pos.z = j["MassComponents"][std::to_string(i)]["position"][2u];

        powermodel power;
        power.off = j["MassComponents"][std::to_string(i)]["powermodel"][0u];
        power.idle = j["MassComponents"][std::to_string(i)]["powermodel"][1u];
        power.use = j["MassComponents"][std::to_string(i)]["powermodel"][2u];
        power.max = j["MassComponents"][std::to_string(i)]["powermodel"][3u];

        quaternion rot;
        component massComponent(mass, moi, pos, rot, power);
        std::string key = j["MassComponents"][std::to_string(i)]["key"];
        returnConfig.components.insert(std::pair<std::string, component>(key, massComponent));
    }

    // parse and map thruster components
    for(int i = 0; i < j["ThrusterComponentNo"]; i++)
    {
        double maxT = j["ThrusterComponents"][std::to_string(i)]["maxThrust"];
        double mTFrac = j["ThrusterComponents"][std::to_string(i)]["minThrustFraction"];
        double tFrac = j["ThrusterComponents"][std::to_string(i)]["thrustFraction"];
        double spImpulse = j["ThrusterComponents"][std::to_string(i)]["specificImpulse"];
        
        bool attitudeControl;
        if (j["ThrusterComponents"][std::to_string(i)]["attitudeControl"] == 0) {
            attitudeControl = false;
        } else {
            attitudeControl = true;
        }
        
        std::array<double, 3> tAxis;
        tAxis[0] = j["ThrusterComponents"][std::to_string(i)]["thrustAxis"][0u];
        tAxis[1] = j["ThrusterComponents"][std::to_string(i)]["thrustAxis"][1u];
        tAxis[2] = j["ThrusterComponents"][std::to_string(i)]["thrustAxis"][2u];
     
        position tCenter;
        tCenter.x = j["ThrusterComponents"][std::to_string(i)]["thrustCenter"][0u];
        tCenter.y = j["ThrusterComponents"][std::to_string(i)]["thrustCenter"][1u];
        tCenter.z = j["ThrusterComponents"][std::to_string(i)]["thrustCenter"][2u];

        thruster ThrusterComponent(maxT, mTFrac, tFrac, spImpulse, attitudeControl, tAxis, tCenter);
        
        std::string key = j["ThrusterComponents"][std::to_string(i)]["key"];
        returnConfig.thrusters.insert (std::pair<std::string, thruster>(key, ThrusterComponent));
        
    }

    // parse and map fuel tanks
    for(int i = 0; i < j["FuelTankNo"]; i++)
    {
        double fuelmass = j["FuelTanks"][std::to_string(i)]["fuelMass"];
        double fuelcap = j["FuelTanks"][std::to_string(i)]["fuelCapacity"];
        std::string fueltype = j["FuelTanks"][std::to_string(i)]["fuelType"];
        position pos;
        pos.x = j["FuelTanks"][std::to_string(i)]["position"][0u];
        pos.y = j["FuelTanks"][std::to_string(i)]["position"][1u];
        pos.z = j["FuelTanks"][std::to_string(i)]["position"][2u];
        fueltank FuelTank(fueltype, fuelmass, fuelcap, pos);

        std::string key = j["FuelTanks"][std::to_string(i)]["key"];

        returnConfig.fueltanks.insert (std::pair<std::string, fueltank>(key, FuelTank));
    }

    // // parse and map rotators
    for(int i = 0; i < j["RotatorNo"]; i++)
    {
        double mDipole = j["Rotators"][std::to_string(i)]["maxDipoleMoment"];
        double torque = j["Rotators"][std::to_string(i)]["torque"];
        double storedMomentum = j["Rotators"][std::to_string(i)]["storedMomentum"];
        double rotSpeed = j["Rotators"][std::to_string(i)]["rotationSpeed"];
        momentofinertia moi;
        moi.Ixx = j["Rotators"][std::to_string(i)]["moi"][0u];
        moi.Iyy = j["Rotators"][std::to_string(i)]["moi"][1u];
        moi.Izz = j["Rotators"][std::to_string(i)]["moi"][2u];
        std::array<double, 3> direction;
        direction[0] = j["Rotators"][std::to_string(i)]["normalVector"][0u];
        direction[1] = j["Rotators"][std::to_string(i)]["normalVector"][1u];
        direction[2] = j["Rotators"][std::to_string(i)]["normalVector"][2u];

        rotator RotatorComponent(mDipole, torque, storedMomentum, rotSpeed, moi, direction);

        std::string key = j["Rotators"][std::to_string(i)]["key"];

        returnConfig.rotators.insert (std::pair<std::string, rotator>(key, RotatorComponent));
    }

    return returnConfig;

};
}

#endif // JSONPARSER_H.begin();