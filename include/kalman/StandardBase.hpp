// The MIT License (MIT)
//
// Copyright (c) 2015-2017 Markus Herb
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
#ifndef KALMAN_STANDARDBASE_HPP_
#define KALMAN_STANDARDBASE_HPP_

#include "Types.hpp"

namespace Kalman {
    
    /**
     * @brief Abstract base class for standard (non-square root) filters and models
     * 
     * @param StateType The vector-type of the system state
     */
    template<class StateType>
    class StandardBase
    {
    public:
        //! Numeric scalar type
        typedef typename StateType::Scalar T;
        //! State Covariance type aliases
        typedef Kalman::Covariance<T, StateType> Covariance;
        //! State Covariance Square Root type aliases
        typedef Kalman::CovarianceSquareRoot<T, StateType> CovarianceSquareRoot;

    protected:
        //! Covariance
        Covariance P;

    public:
        /**
         * @brief Get covariance
         */
        const Covariance& getCovariance() const
        {
            return P;
        }

        /**
         * @brief Get covariance (as square root)
         */
        CovarianceSquareRoot getCovarianceSquareRoot() const
        {
            return CovarianceSquareRoot(P);
        }

        /**
         * @brief Set Covariance
         */
        bool setCovariance(const Covariance& covariance)
        {
            P = covariance;
            return true;
        }

    protected:
        StandardBase()
        {
            P.setIdentity();
        }
    };
}

#endif
