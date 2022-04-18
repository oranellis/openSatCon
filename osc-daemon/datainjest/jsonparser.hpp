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

craftconfig parseJson(std::string jsonPath) {
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
        moi.Ixx = j["MassComponents"][std::to_string(i)]["moi"][0u];
        moi.Iyy = j["MassComponents"][std::to_string(i)]["moi"][1u];
        moi.Izz = j["MassComponents"][std::to_string(i)]["moi"][2u];
        position pos;
        pos.x = j["MassComponents"][std::to_string(i)]["position"][0u];
        pos.y = j["MassComponents"][std::to_string(i)]["position"][1u];
        pos.z = j["MassComponents"][std::to_string(i)]["position"][2u];
        powermodel power;
        std::cout << j["MassComponents"][std::to_string(i)]["powermodel"];
        quaternion rot;
        component massComponent(mass, moi, pos, rot, power);

        returnConfig.components.insert(std::pair<std::string, component>(j["MassComponents"][std::to_string(i)]["key"], massComponent));
    }

    // parse and map thruster components
    // for(int i = 0; i < j["ThrustComponentNo"]; i++)
    // {
    //     thruster thrustComponent;
    //     thrustComponent.maxThrust = j["MassComponents"][std::to_string(i)]["maxThrust"];
    //     thrustComponent.minThrustFraction = j["MassComponents"][std::to_string(i)]["minThrustFraction"];
    //     thrustComponent.thrustFraction = j["MassComponents"][std::to_string(i)]["ThrustFraction"];
    //     thrustComponent.specificImpulse = j["MassComponents"][std::to_string(i)]["specificImpulse"];
    //     thrustComponent.thrustAxis = j["MassComponents"][std::to_string(i)]["thrustAxis"];

    //     returnConfig.thrusters.insert (std::pair<std::string, thruster>(j["ThrustComponents"][std::to_string(i)]["key"], thrustComponent));
        
    // }

    // // parse and map fuel tanks
    // for(int i = 0; i < j["FuelTankNo"]; i++)
    // {
    //     fueltank FuelTank;
    //     FuelTank.fuelMass = j["FuelTanks"][std::to_string(i)]["FuelMass"];
    //     FuelTank.fuelMassPos = j["FuelTanks"][std::to_string(i)]["position"];

    //     returnConfig.fueltanks.insert (std::pair<std::string, fueltank>(j["FuelTanks"][std::to_string(i)]["key"], FuelTank));

    //     // returnConfig.fueltanks.insert (std::pair<std::string, double>(j["FuelTanks"][std::to_string(i)]["key"], j["FuelTanks"][std::to_string(i)]["FuelMass"]));
    //     // returnConfig.fueltanks.insert (std::pair<std::string, position>(j["FuelTanks"][std::to_string(i)]["key"], j["FuelTanks"][std::to_string(i)]["position"]));
    // }

    // // parse and map rotators
    // for(int i = 0; i < j["ReactionWheelNo"]; i++)
    // {
    //     returnConfig.rotators.insert (std::pair<std::string, double>(j["Rotators"][std::to_string(i)]["key"], j["Rotators"][std::to_string(i)]["Torque"]));
    //     returnConfig.rotators.insert (std::pair<std::string, double>(j["Rotators"][std::to_string(i)]["key"], j["Rotators"][std::to_string(i)]["RotationSpeed"]));
    // }

    return returnConfig;
// access data via craft["ThrusterComponents"]["normal"] etc
};
}

#endif // JSONPARSER_H