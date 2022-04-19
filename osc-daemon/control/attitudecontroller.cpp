#include "attitudecontroller.hpp"

namespace osc::attitudecontrol {

    const int CONTROL_LOOP_FREQ = 8000;
    const double KP = 0.01;

    void controlLoopThread(bool *interupt, task *curTask, craftcontroller *controller) {
        vec3 allignmentAxis = vec3(1, 0, 0);


        while (!*interupt) {

            auto loopStartTime = std::chrono::steady_clock::now();
            vec3 setpointAxis = curTask->getPointingDirection();
            vec3 currentAxis = (controller->getRotation()).rotate(allignmentAxis);

            std::array<double, 4> axisAngleSetpoint = quaternion(currentAxis, setpointAxis).getAxisAngle();

            vec3 scaledAngularAccelSetpoint = { axisAngleSetpoint[0] * axisAngleSetpoint[1], axisAngleSetpoint[0] * axisAngleSetpoint[2], axisAngleSetpoint[0] * axisAngleSetpoint[3] };
            vec3 angularAccelSetpoint = scaledAngularAccelSetpoint/2*M_PI;

            vec3 momentCommand = angularAccelSetpoint * KP; // Proportional controller, additional terms required for full PID
            ftModel ftCommand = {0, 0, 0, momentCommand[0], momentCommand[1], momentCommand[2]};

            controller->forcesToCommands(ftCommand);

            std::chrono::time_point suspendUntil = loopStartTime + std::chrono::microseconds(1000000/CONTROL_LOOP_FREQ); // Period represented in microseconds, sets the time to start the next control loop
            std::this_thread::sleep_until(loopStartTime); // Will run slow if loop takes longer than 250us
        }
    }
}