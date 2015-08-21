cmake_minimum_required(VERSION 2.8)
project(Kalman)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(Eigen3 REQUIRED)

set (HEADERS
    # Base
    "include/kalman/KalmanFilterBase.h"
    "include/kalman/StandardBase.h"
    "include/kalman/SquareRootBase.h"
    "include/kalman/StandardFilterBase.h"
    "include/kalman/SquareRootFilterBase.h"
    
    # EKF
    "include/kalman/ExtendedKalmanFilter.h"
    "include/kalman/SquareRootExtendedKalmanFilter.h"
    
    # UKF
    "include/kalman/UnscentedKalmanFilterBase.h"
    "include/kalman/UnscentedKalmanFilter.h"
    "include/kalman/SquareRootUnscentedKalmanFilter.h"
    
    # Models
    "include/kalman/MeasurementModel.h"
    "include/kalman/SystemModel.h"
    "include/kalman/LinearizedMeasurementModel.h"
    "include/kalman/LinearizedSystemModel.h"
    
    # Utils
    "include/kalman/Matrix.h"
    "include/kalman/Types.h"
)

include_directories(
    "include"
    ${Eigen3_INCLUDE_DIR}
)

# Flags
set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -pedantic -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DEIGEN_NO_DEBUG")

# Examples
add_executable(example_robot1 "examples/Robot1/main.cpp")

