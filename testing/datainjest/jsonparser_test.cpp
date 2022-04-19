#include "../../osc-daemon/datainjest/jsonparser.hpp"
#include <iostream>
#include <fstream>

using json = nlohmann::json;

int main() {

    std::string path = "../../osc-daemon/datainjest/examplecraft.json";
  
    osc::craftconfig konfig; 

    konfig = osc::parseJson(path);

    for ( auto p : konfig.components)
    {
        std::cout << p.first << "\n";
    }
        for ( auto p : konfig.thrusters)
    {
        std::cout << p.first << "\n";
    }
        for ( auto p : konfig.fueltanks)
    {
        std::cout << p.first << "\n";
    }
        for ( auto p : konfig.rotators)
    {
        std::cout << p.first << "\n";
    }
    return 0;
}