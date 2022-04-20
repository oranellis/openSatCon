#ifndef SCHEDULER
#define SCHEDULER

#include "task.hpp"

#include <vector>

namespace osc {

    /** \class scheduler
    \brief handles task processing
    scheduler class to define an object to handle task processing
    */
    class scheduler {

        private:

        /// @param state state of the scheduler, true being active
        bool state = true;
        /// @param ptasks vector of tasks 
        std::vector<task> * ptasks = new std::vector<task>;

        public:

        /** \fn active()
        returns the state of the scheduler 
        */
        bool active() {
            return state;
        }

        /** getNext() returns the index of the next task 
        @param[out] closestTaskIndex output index of the next task
        */
        task getNext() {
            
            std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
            int closestTaskIndex = -1;

            if (ptasks->size()>0) {

                std::chrono::time_point<std::chrono::system_clock> closestTime = std::chrono::time_point<std::chrono::system_clock>::max();

                for (int i = 0; i < ptasks->size(); i++) {
                    if ((ptasks->at(i)).getStartTime() > now && ptasks->at(i).getStartTime() < closestTime) {
                        closestTime = ptasks->at(i).getStartTime();
                        closestTaskIndex = i;
                    }
                }
                if (closestTaskIndex >= 0) {
                    return ptasks->at(closestTaskIndex);
                }
                return task();
            }
        }

        /** \fn addTask(argTask) adds \p argTask into the tasks list
        @param[in] argTask task to add
        
        */
        void addTask(task argTask) {
            ptasks->push_back(argTask);
        }

    };
} // namespace osc

#endif