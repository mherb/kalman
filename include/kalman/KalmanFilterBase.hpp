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
#ifndef KALMAN_KALMANFILTERBASE_HPP_
#define KALMAN_KALMANFILTERBASE_HPP_

#include "Matrix.hpp"
#include "Types.hpp"

namespace Kalman {
    
    /**
     * @brief Abstract base class for all Kalman Filters
     * 
     * @param StateType The vector-type of the system state (usually some type derived from Kalman::Vector)
     */
    template<class StateType>
    class KalmanFilterBase
    {
    public:
        static_assert(/*StateType::RowsAtCompileTime == Dynamic ||*/StateType::RowsAtCompileTime > 0,
                      "State vector must contain at least 1 element" /* or be dynamic */);
        static_assert(StateType::ColsAtCompileTime == 1, "State type must be a column vector");

        //! Numeric scalar type
        typedef typename StateType::Scalar T;
        
        //! Type of the state vector
        typedef StateType State;
        
    protected:
        //! Estimated state
        State x;
        
    public:
        /**
         * Get current state estimate
         */
        const State& getState() const
        {
            return x;
        }
        
        /**
         * @brief Initialize state
         * @param initialState The initial state of the system
         */
        void init(const State& initialState)
        {
            x = initialState;
        }
    protected:
        /**
         * @brief Protected constructor
         */
        KalmanFilterBase()
        {
        }
    };
}

#endif
