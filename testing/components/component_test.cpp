#include "../../osc-daemon/components/component.hpp"
#include "../../osc-daemon/osctypes.hpp"
#include <iostream>

int main() {
    osc::momentofinertia moi;
    moi.Ixx = 1; moi.Iyy = 3; moi.Izz = 0.2;
    osc::component comp(12,moi,osc::position(0.2,0.1,0.3), osc::quaternion(osc::vec3(1,0,0), osc::vec3(0,1,0)), osc::powermodel());
    std::cout << comp.getMass() << std::endl;
    return 0;
}