#ifndef SCHEDULER
#define SCHEDULER

#include "task.hpp"

#include <vector>

namespace osc {
    class scheduler {

        private:

        bool state = true;
        std::vector<task> * ptasks = new std::vector<task>;

        public:

        bool active() {
            return state;
        }

        int getNext() {
            
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
                
                return closestTaskIndex;
            }
        }

        void addTask(task argTask) {
            ptasks->push_back(argTask);
        }

    };
} // namespace osc

#endif