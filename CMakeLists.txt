cmake_minimum_required(VERSION 3.22)

project("osc-daemon") 

set(SOURCE_FILES src/main.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})