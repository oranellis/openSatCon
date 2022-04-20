#include "craftcontroller.hpp"

namespace osc {

    /** \fn initModel()
    Initialises the model, returning true when the model initialises 
    */
    bool craftcontroller::initModel() {
        std::string pathString;
        // std::cout << "Enter path to craft configuration: ";
        // std::cin >> pathString;

        std::cout << "Using the default craft config " << std::endl;
        pathString = "../osc-daemon/datainjest/examemecraft.json";

        craftconfig config = parseJson(pathString);
        
        if (!config.populated()) return false;
        return true;
    }

    /** /fn getMaxThrust()
    returns the \p maxThrust parameter
    */
    double craftcontroller::getMaxThrust() {
        return maxThrust;
    }

    /** /fn getTransferISP()
    returns the \p transferISP parameter
    */
    double craftcontroller::getTransferISP() {
        return transferISP;
    }

    /** \fn recomputeComponentDeps()
    recomputes craft parameters from a change in component configuration
    */
    void craftcontroller::recomputeComponentDeps() {

        cg = position(0,0,0);
        moi = momentofinertia();
        mass = 0;
        wetMass = 0;

        for (auto i=components.begin(); i!=components.end(); i++) {
            cg = ((i->second.getPos().multiply(i->second.getMass())) + \
                    (cg.multiply((mass + wetMass)))).divide \
                    (i->second.getMass() + mass + wetMass);
            mass += i->second.getMass();
        }
        
        for (auto i=components.begin(); i!=components.end(); i++) {
            component* thisComponent = &i->second; // Assign pointer to class for readability
            position offset = thisComponent->getPos() - cg;
            moi.addMass(thisComponent->getMoi(),thisComponent->getMass(),offset);
        }

        for (auto i=fueltanks.begin(); i!=fueltanks.end(); i++) {
            cg = ((i->second.getFuelPos().multiply(i->second.getFuelMass())) + \
                    (cg.multiply((mass + wetMass)))).divide \
                    (i->second.getFuelMass() + mass + wetMass);
            wetMass += i->second.getFuelMass();
        }

        for (auto i=thrusters.begin(); i!=thrusters.end(); i++) {
            if (i->second.getAttitudeControl()) {
                attitudeActuators.insert({i->first, ftModel(i->second.getThrustAxis()*i->second.getMaxThrust(), components[i->first].getPos())}); // Add ftModel and string identifier to attitudeActuators
            }
        }

        for (auto i=rotators.begin(); i!=rotators.end(); i++) {
            vec3 torqueVec = i->second.getDirection()*i->second.getTorque();
            attitudeActuators.insert({i->first, ftModel(0, 0, 0, torqueVec[0], torqueVec[1], torqueVec[2])}); // Add ftModel and string identifier to attitudeActuators
            std::cout << "_Debug_ rotator torque in x: " << torqueVec[0] << std::endl;
        }
    }

    /** \fn forcesToCommands()
    @param[in] setpoint input an ftModel of the setpoint
    returns a set of thruster commands
    */
    void craftcontroller::forcesToCommands(ftModel setpoint) {

        std::map<std::string, double> currThrusterCommands;

        bool matched = false;
        ftModel sp = setpoint;
        int count = 0;

        while (!matched) {
            int currAxis = sp.getDominantAxis();
            std::string actuatorID;
            for (auto i=attitudeActuators.begin(); i!= attitudeActuators.end(); i++) {
                if (i->second.getDominantAxis() == currAxis) {
                    double command = sp[currAxis]/i->second[currAxis]; // Eliminate command error in dominant axis

                    if (currThrusterCommands.count(i->first)>0) {
                        currThrusterCommands[i->first] = currThrusterCommands[i->first] + command;
                    }

                    else {
                        currThrusterCommands.insert({i->first, command});
                    }
                }
                if (sp.mag()<0.01*setpoint.mag() || count > 100) {
                    matched = true;
                }
            }

            thrusterCommands = currThrusterCommands;
        }
    }

    /** \fn controlLoopThread()
    the function that is run in a thread responsible for generating the actuator outputs based on the current task
    */
    void craftcontroller::controlLoopThread() {

        vec3 allignmentAxis = vec3(1, 0, 0);

        while (!taskInterupt) {

            auto loopStartTime = std::chrono::steady_clock::now();
            vec3 setpointAxis = currTask.getPointingDirection();
            vec3 currentAxis = rotation.rotate(allignmentAxis);

            std::array<double, 4> axisAngleSetpoint = quaternion(currentAxis, setpointAxis).getAxisAngle();

            vec3 scaledAngularAccelSetpoint = { axisAngleSetpoint[0] * axisAngleSetpoint[1], axisAngleSetpoint[0] * axisAngleSetpoint[2], axisAngleSetpoint[0] * axisAngleSetpoint[3] };
            vec3 angularAccelSetpoint = scaledAngularAccelSetpoint/2*M_PI;

            vec3 momentCommand = angularAccelSetpoint * KP; // Proportional controller, additional terms required for full PID
            ftModel ftCommand = {0, 0, 0, momentCommand[0], momentCommand[1], momentCommand[2]};

            forcesToCommands(ftCommand);

            auto suspendUntil = loopStartTime + std::chrono::microseconds(1000000/CONTROL_LOOP_FREQ); // Period represented in microseconds, sets the time to start the next control loop
            std::this_thread::sleep_until(loopStartTime); // Will run slow if loop takes longer than 125us
        }
    }

    void craftcontroller::outputThread() {

        wiringPiSetup();

        pinMode(PIN_X_POS, OUTPUT);
        pinMode(PIN_X_NEG, OUTPUT);
        pinMode(PIN_Y_POS, OUTPUT);
        pinMode(PIN_Y_NEG, OUTPUT);
        pinMode(PIN_Z_POS, OUTPUT);
        pinMode(PIN_Z_NEG, OUTPUT);
        pinMode(PIN_T, OUTPUT);

        while (!taskInterupt) {

            auto timeNow = std::chrono::steady_clock::now();
            digitalWrite(PIN_X_POS, 1);
            digitalWrite(PIN_X_NEG, 1);
            digitalWrite(PIN_Y_POS, 1);
            digitalWrite(PIN_Y_NEG, 1);
            digitalWrite(PIN_Z_POS, 1);
            digitalWrite(PIN_Z_NEG, 1);
            digitalWrite(PIN_T, 1);
            std::this_thread::sleep_until(timeNow+std::chrono::milliseconds(500));

            digitalWrite(PIN_X_POS, 0);
            digitalWrite(PIN_X_NEG, 0);
            digitalWrite(PIN_Y_POS, 0);
            digitalWrite(PIN_Y_NEG, 0);
            digitalWrite(PIN_Z_POS, 0);
            digitalWrite(PIN_Z_NEG, 0);
            digitalWrite(PIN_T, 0);
            std::this_thread::sleep_until(timeNow+std::chrono::milliseconds(1000));

        }
        
    }

    /// If model fails to initialise
    craftcontroller::craftcontroller() {
        if (!initModel()) throw std::runtime_error("Craft model failed to initialise");
    }

    /** \fn beginControl()
    control begins and scheduler is activated 
    */
    void craftcontroller::beginExampleControl() {
        // scheduler schedule = scheduler();
        std::cout << "Beginning control" << std::endl;
        
        // while(schedule.active()) {
        //     task currTask = schedule.getNext();
            
        //     /// @param interupt Acts as a control flag to stop the attitude control thread
        //     bool interupt = false;

        //     std::thread attitudeControlThread(controlLoopThread, &currTask, &interupt);
        // }

        currTask = task(taskType::example);

        std::cout << "Example task generated" << std::endl;

        std::thread outputThread(&craftcontroller::controlLoopThread, this);

        // std::thread attitudeControlThread(&craftcontroller::controlLoopThread, this);
    }
};