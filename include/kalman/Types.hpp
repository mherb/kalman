// The MIT License (MIT)
//
// Copyright (c) 2015 Markus Herb
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
#ifndef KALMAN_TYPES_HPP_
#define KALMAN_TYPES_HPP_

#include "Matrix.hpp"

namespace Kalman
{
    /**
     * @class Kalman::SquareMatrix
     * @brief Template type representing a square matrix
     * @param T The numeric scalar type
     * @param N The dimensionality of the Matrix
     */
    template<typename T, int N>
    using SquareMatrix = Matrix<T, N, N>;
    
    /**
     * @class Kalman::Covariance
     * @brief Template type for covariance matrices
     * @param Type The vector type for which to generate a covariance (usually a state or measurement type)
     */
    template<class Type>
    using Covariance = SquareMatrix<typename Type::Scalar, Type::RowsAtCompileTime>;
    
    /**
     * @class Kalman::CovarianceSquareRoot
     * @brief Template type for covariance square roots
     * @param Type The vector type for which to generate a covariance (usually a state or measurement type)
     */
    template<class Type>
    using CovarianceSquareRoot = Cholesky< Covariance<Type> >;
    
    /**
     * @class Kalman::KalmanGain
     * @brief Template type of Kalman Gain
     * @param State The system state type
     * @param Measurement The measurement type
     */
    template<class State, class Measurement>
    using KalmanGain = Matrix<typename State::Scalar,
                              State::RowsAtCompileTime,
                              Measurement::RowsAtCompileTime>;
    
    /**
     * @class Kalman::Jacobian
     * @brief Template type of jacobian of A w.r.t. B
     */
    template<class A, class B>
    using Jacobian = Matrix<typename A::Scalar,
                            A::RowsAtCompileTime,
                            B::RowsAtCompileTime>;
}

#endif
