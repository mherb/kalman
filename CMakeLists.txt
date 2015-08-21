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

# Unit Tests
set(TESTS
    # EKF
    "test/ExtendedKalmanFilter.cpp"
    
    # UKF
    "test/UnscentedKalmanFilterBase.cpp"
    "test/UnscentedKalmanFilter.cpp"
    "test/SquareRootUnscentedKalmanFilter.cpp"
    
    # Utils
    "test/Matrix.cpp"
)

enable_testing()
add_subdirectory("external/googletest")
include_directories("external/googletest/include")
add_executable(kalman_test ${TESTS})
target_link_libraries(kalman_test PRIVATE gtest gtest_main)
add_test(Kalman kalman_test)
