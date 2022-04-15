#include <iostream>
#include "../../includes/json.hpp"
#include <list>

using namespace std;
using json = nlohmann::json;

// change name to parse()
int main() 
{
    std::cout << "start";
    map<string, std::list<double> > Components;
    
    json j = json::parse("craft.json");
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

    std::cout << "Components";
    //parse thrust components
}
// take in json data 
// take a list of components each with component.cpp's data types
// map via ID name

// access data via craft["ThrusterComponents"]["normal"] etc
