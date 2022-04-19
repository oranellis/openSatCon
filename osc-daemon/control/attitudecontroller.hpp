#ifndef ATTITUDECONTROLLER_H
#define ATTITUDECONTROLLER_H

#include <chrono>

#include "../osctypes.hpp"
#include "task.hpp"

namespace osc::attitudecontrol {

    void controlLoopThread(bool *interupt, task *curTask) {}

}

#endif // ATTITUDECONTROLLER_H