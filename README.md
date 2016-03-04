# Kalman Filter Library

[![Build Status](https://travis-ci.org/mherb/kalman.svg?branch=master)](https://travis-ci.org/mherb/kalman)

This is a header-only C++11 library implementing common variants of the well-known [Kalman-Filter](https://en.wikipedia.org/wiki/Kalman_filter).
Currently implementations of these filter variants are included:

* Extended Kalman Filter (EKF)
* Square Root Extended Kalman Filter (SR-EKF)
* Unscented Kalman Filter (UKF)
* Square Root Unscented Kalman Filter (SR-UKF)

## Dependencies

This library makes heavy use of the excellent [Eigen3 library](http://eigen.tuxfamily.org) for linear algebra operations and is thus a required dependency.

## Usage
In order to use the library to do state estimation, a number of things have to be done:

1. Define a state-vector type
2. (Optional) Define a control-vector type
3. Define a system model
4. Define one (or more) measurement models with corresponding measurement vector types

### Example
A fairly worked out example on how to use the library is given in `examples/Robot1` with detailed commentary.

### State Vector
The state vector defines the state variables of your system that should be estimated.
You can use the readily available `Kalman::Vector` template type as your vector or derive your own specialized state vector from that.

### Control Vector
In case your system has some control input, a control vector has to be defined analogously to the state vector.

### System Model
The system model defines how the system state evolves over time, i.e. from one time-step to the next given some control input.
The transition function is in general non-linear. Any system model must derive from the base `SystemModel` class template.
In case a linearized filter such as the Extended Kalman Filter should be used, then the system model must be given as linearized model by deriving from `LinearizedSystemModel` and defining the corresponding jacobians.

Note that linearized models can of course also be used with fully non-linear filters such as the Unscented Kalman Filter.

### Measurement Vector
The measurement vector represents the measurement taken by some sensors and has to be defined analogously to the state and control vectors.

### Measurement Model
The measurement model defines how a measurement is related to the system state, i.e. it maps a system state to the expected sensor measurement.
Measurement models must derive from the class template `MeasurementModel` or, in case of linearized models for EKFs, from `LinearizedMeasurementModel`.

## FAQ

__The filters are running very slowly, why is that and how can I make them faster?__
By default, operations in Eigen include a lot of debug code, such as checking for valid matrix and vector bounds and so on.
To speed things up, these checks can be disabled using the pre-processor define

    -DEIGEN_NO_DEBUG

which is also automatically set when using the general

    -DNDEBUG

flag. In addition to that the regular optimization flags including `-O2` will make things faster.


## License

The MIT License (MIT)

Copyright (c) 2015 [mherb](https://github.com/mherb)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
