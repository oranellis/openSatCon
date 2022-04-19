#include "attitudecontroller.hpp"

namespace osc::attitudecontrol {

    void controlLoopThread(bool *interupt, task *curTask) {
        vec3 allignmentAxis = vec3(1, 0, 0);
        vec3 setpointAxis = curTask->getPointingDirection();
    }

}