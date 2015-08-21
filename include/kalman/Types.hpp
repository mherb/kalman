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
    using Covariance = SquareMatrix<typename Type::Scalar, Type::length>;
    
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
    using KalmanGain = Matrix<typename State::Scalar, State::length, Measurement::length>;
    
    /**
     * @class Kalman::Jacobian
     * @brief Template type of jacobian of A w.r.t. B
     */
    template<class A, class B>
    using Jacobian = Matrix<typename A::Scalar, A::length, B::length>;
}

#endif
