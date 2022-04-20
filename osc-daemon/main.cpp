#include <iostream>

#include "control/craftcontroller.hpp"

int main(int argc, char** args) {

    std::cout << "Welcome to Open Satelite Control" << std::endl;

    osc::craftcontroller controller;

    controller.beginExampleControl();

    std::cout << "Finished execution, closing." << std::endl;

    return 0;

}