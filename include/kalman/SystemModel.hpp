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
#ifndef KALMAN_SYSTEMMODEL_HPP_
#define KALMAN_SYSTEMMODEL_HPP_

#include <type_traits>

#include "Matrix.hpp"
#include "StandardBase.hpp"

namespace Kalman {
    /**
     * @brief Abstract base class of all system models
     *
     * @param StateType The vector-type of the system state (usually some type derived from Kalman::Vector)
     * @param ControlType The vector-type of the control input (usually some type derived from Kalman::Vector)
     * @param CovarianceBase The class template used for covariance storage (must be either StandardBase or SquareRootBase)
     */
    template<class StateType, class ControlType = Vector<typename StateType::Scalar, 0>, template<class> class CovarianceBase = StandardBase>
    class SystemModel : public CovarianceBase<StateType>
    {
        static_assert(/*StateType::RowsAtCompileTime == Dynamic ||*/ StateType::RowsAtCompileTime > 0,
                      "State vector must contain at least 1 element" /* or be dynamic */);
        static_assert(/*ControlType::RowsAtCompileTime == Dynamic ||*/ ControlType::RowsAtCompileTime >= 0,
                      "Control vector must contain at least 0 elements" /* or be dynamic */);
        static_assert(std::is_same<typename StateType::Scalar, typename ControlType::Scalar>::value,
                      "State and Control scalar types must be identical");
    public:
        //! System state type
        typedef StateType State;
        
        //! System control input type
        typedef ControlType Control;
        
    public:
        /**
         * @brief State Transition Function f
         * 
         * Computes the predicted system state in the next timestep given
         * the current state x and the control input u
         */
        virtual State f(const State& x, const Control& u) const = 0;
        
    protected:
        SystemModel() {}
        virtual ~SystemModel() {}
    };
}

#endif
