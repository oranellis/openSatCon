cmake_minimum_required(VERSION 3.22)

project("osc-daemon") 

set(SOURCE_FILES    osc-daemon/main.cpp
                    osc-daemon/control/craftcontroller.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})