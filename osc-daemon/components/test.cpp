#include "component.hpp"

int main () {
    osc::component component(12.3, osc::position(0,0,0), osc::powermodel());
    component.testPrint();
}