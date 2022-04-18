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
        returnConfig.components.insert (std::pair<std::string, component.mass>(j["MassComponents"][std::to_string(i)]["key"], j["MassComponents"][std::to_string(i)]["mass"]));
        returnConfig.components.insert (std::pair<std::string, component.moi>(j["MassComponents"][std::to_string(i)]["key"], j["MassComponents"][std::to_string(i)]["moi"]));
        returnConfig.components.insert (std::pair<std::string, component.pos>(j["MassComponents"][std::to_string(i)]["key"], j["MassComponents"][std::to_string(i)]["position"]));
        returnConfig.components.insert (std::pair<std::string, component.power>(j["MassComponents"][std::to_string(i)]["key"], j["MassComponents"][std::to_string(i)]["powermodel"]));

        std::cout << returnConfig.components["thruster"];
    }

    // parse and map thruster components
    for(int i = 0; i < j["ThrustComponentNo"]; i++)
    {
        returnConfig.thrusters.insert (std::pair<std::string, thruster.maxThrust>(j["ThrustComponents"][std::to_string(i)]["key"], j["MassComponents"][std::to_string(i)]["maxThrust"]));
        returnConfig.thrusters.insert (std::pair<std::string, thruster.minThrustFraction>(j["ThrustComponents"][std::to_string(i)]["key"], j["MassComponents"][std::to_string(i)]["minThrustFraction"]));
        returnConfig.thrusters.insert (std::pair<std::string, thruster.ThrustFraction>(j["ThrustComponents"][std::to_string(i)]["key"], j["MassComponents"][std::to_string(i)]["ThrustFraction"]));
        returnConfig.thrusters.insert (std::pair<std::string, thruster.specificImpulse>(j["ThrustComponents"][std::to_string(i)]["key"], j["MassComponents"][std::to_string(i)]["specificImpulse"]));
        returnConfig.thrusters.insert (std::pair<std::string, thruster.thrustAxis>(j["ThrustComponents"][std::to_string(i)]["key"], j["MassComponents"][std::to_string(i)]["thrustAxis"]));

        std::cout << returnConfig.components["thruster"];
    }

    // parse and map fuel tanks
    for(int i = 0; i < j["FuelTankNo"]; i++)
    {
        returnConfig.fueltanks.insert (std::pair<std::string, fueltank.fuelMass>(j["FuelTanks"][std::to_string(i)]["key"], j["FuelTanks"][std::to_string(i)]["FuelMass"]));
        returnConfig.fueltanks.insert (std::pair<std::string, fueltank.fuelMassPos>(j["FuelTanks"][std::to_string(i)]["key"], j["FuelTanks"][std::to_string(i)]["position"]));

        std::cout << returnConfig.components["fuel_tank"];
    }

    // parse and map rotators
    for(int i = 0; i < j["ReactionWheelNo"]; i++)
    {
        returnConfig.rotators.insert (std::pair<std::string, rotator.reactionWheel.Torque>(j["Rotators"][std::to_string(i)]["key"], j["Rotators"][std::to_string(i)]["Torque"]));
        returnConfig.rotators.insert (std::pair<std::string, rotator.reactionWheel.RotationSpeed>(j["Rotators"][std::to_string(i)]["key"], j["Rotators"][std::to_string(i)]["RotationSpeed"]));

        std::cout << returnConfig.components["rotator_1"];
    }

// access data via craft["ThrusterComponents"]["normal"] etc
}

}

#endif // JSONPARSER_H