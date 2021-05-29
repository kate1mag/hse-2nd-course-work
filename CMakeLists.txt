cmake_minimum_required(VERSION 3.15)
project(panther_master3)

set(CMAKE_CXX_STANDARD 17)

add_executable(panther_master CoordinateDescentMethod/CoordDescMethod_test.cpp
        CoordinateDescentMethod/CoordinateDescentMethod.h CoordinateDescentMethod/parallCoordDescMethod.h
        CoordinateDescentMethod/CoordinateDescentMethod_old.h)
include_directories(.)
add_subdirectory(Google_tests)