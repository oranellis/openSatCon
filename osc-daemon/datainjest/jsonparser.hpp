#ifndef JSONPARSER_H
#define JSONPARSER_H

#include <iostream>
#include "../../includes/json.hpp"
#include "../components/component.hpp"
#include "../components/fueltank.hpp"
#include "../components/actuators/thruster.hpp"
#include "../components/actuators/rotator.hpp"
#include <list>
#include <vector>
#include <fstream>
//#include "../components/component.hpp"

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
    map<string, std::list<double> > Components;
    
    json j = json::parse(jsonPath);
    // std::ifstream ifs("craft.json");
    // json j = json::parse(ifs);

    // std::string str(R"({"json": "beta"})");
    // json js = json::parse(str);

    //parse mass components
    for(int i = 0; i < 2; i++)
    {
        //access data through j["name"]["name"]
        Components.insert(std::make_pair(std::to_string(i), j["MassComponents"][std::to_string(i)]["componentType"]));
        Components.insert(std::make_pair(std::to_string(i), j["MassComponents"][std::to_string(i)]["mass"]));
        Components.insert(std::make_pair(std::to_string(i), j["MassComponents"][std::to_string(i)]["moi"]));
        Components.insert(std::make_pair(std::to_string(i), j["MassComponents"][std::to_string(i)]["position"]));
        Components.insert(std::make_pair(std::to_string(i), j["MassComponents"][std::to_string(i)]["powermodel"]));
    }

    // std::cout << "Components"; // debug
    //parse thrust components
}
// take in json data 
// take a list of components each with component.cpp's data types
// map via ID name

// access data via craft["ThrusterComponents"]["normal"] etc
}

#endif // JSONPARSER_H