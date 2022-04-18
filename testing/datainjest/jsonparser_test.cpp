#include "../../osc-daemon/datainjest/jsonparser.hpp"
#include "iostream"

int main() {
    osc::parseJson("../../osc-daemon/datainjest/examplecraft.json");
    return 0;
    //std::cout << returnConfig.component["plate_3"]
}