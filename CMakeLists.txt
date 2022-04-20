cmake_minimum_required(VERSION 3.18)

project("osc-daemon") 

find_package(WiringPi REQUIRED)
find_package(Threads REQUIRED) 

include_directories(${WIRINGPI_INCLUDE_DIRS}) 

set(SOURCE_FILES    osc-daemon/main.cpp
                    osc-daemon/control/craftcontroller.cpp)

add_compile_options(-lwiringPi -lpthread)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME} ${WIRINGPI_LIBRARIES}) 
target_link_libraries(${PROJECT_NAME} ${CMAKE_THREAD_LIBS_INIT})