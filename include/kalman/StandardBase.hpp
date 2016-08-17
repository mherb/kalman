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
#ifndef KALMAN_STANDARDBASE_HPP_
#define KALMAN_STANDARDBASE_HPP_

#include "Types.hpp"

namespace Kalman {
    
    /**
     * @brief Abstract base class for standard (non-square root) filters and models
     * 
     * @param StateType The vector-type of the system state (usually some type derived from Kalman::Vector)
     */
    template<class StateType>
    class StandardBase
    {
    protected:
        //! Covariance
        Covariance<StateType> P;
        
    public:
        /**
         * Get covariance
         */
        const Covariance<StateType>& getCovariance() const
        {
            return P;
        }
        
        /**
         * Get covariance (as square root)
         */
        CovarianceSquareRoot<StateType> getCovarianceSquareRoot() const
        {
            return CovarianceSquareRoot<StateType>(P);
        }
        
        /**
         * Set Covariance
         */
        bool setCovariance(const Covariance<StateType>& covariance)
        {
            P = covariance;
            return true;
        }

        /**
         * @brief Set Covariance using Square Root
         *
         * @param covSquareRoot Lower triangular Matrix representing the covariance
         *                      square root (i.e. P = LLˆT)
         */
        bool setCovarianceSquareRoot(const Covariance<StateType>& covSquareRoot)
        {
            CovarianceSquareRoot<StateType> S;
            S.setL(covSquareRoot);
            P = S.reconstructedMatrix();
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
