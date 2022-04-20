cmake_minimum_required(VERSION 3.18)

project("osc-daemon") 

set(SOURCE_FILES    osc-daemon/main.cpp
                    osc-daemon/control/craftcontroller.cpp)

add_compile_options(-I/usr/local/include -L/usr/local/lib -lwiringPi)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})