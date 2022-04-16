#ifndef TASK
#define TASK

#include <chrono>
#include <iostream>
#include <array>

namespace osc {
    class task {

        private:

        int priority;
        std::chrono::microseconds actionDuration;
        std::array<double, 3>* pointingVec; // this needs to be constantly changing throughout the burn so needs to be dynamic. Maybe make a seperate manoevure object that can be called for the current pointing vector that also includes the endpoint

        public:

    };
}

#endif