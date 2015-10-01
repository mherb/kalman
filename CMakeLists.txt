cmake_minimum_required(VERSION 2.8)
project(Kalman)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(Eigen3 REQUIRED)

set(HEADERS
    # Base
    include/kalman/KalmanFilterBase.hpp
    include/kalman/StandardBase.hpp
    include/kalman/SquareRootBase.hpp
    include/kalman/StandardFilterBase.hpp
    include/kalman/SquareRootFilterBase.hpp

    # EKF
    include/kalman/ExtendedKalmanFilter.hpp
    include/kalman/SquareRootExtendedKalmanFilter.hpp

    # UKF
    include/kalman/UnscentedKalmanFilterBase.hpp
    include/kalman/UnscentedKalmanFilter.hpp
    include/kalman/SquareRootUnscentedKalmanFilter.hpp

    # Models
    include/kalman/MeasurementModel.hpp
    include/kalman/SystemModel.hpp
    include/kalman/LinearizedMeasurementModel.hpp
    include/kalman/LinearizedSystemModel.hpp

    # Utils
    include/kalman/Matrix.hpp
    include/kalman/Types.hpp
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
    test/ExtendedKalmanFilter.cpp
    
    # UKF
    test/UnscentedKalmanFilterBase.cpp
    test/UnscentedKalmanFilter.cpp
    test/SquareRootUnscentedKalmanFilter.cpp
    
    # Utils
    test/Matrix.cpp
)

enable_testing()
add_subdirectory("external/googletest")
include_directories("external/googletest/include")
add_executable(kalman_test ${TESTS} ${HEADERS})
target_link_libraries(kalman_test PRIVATE gtest gtest_main)
add_test(Kalman kalman_test)
