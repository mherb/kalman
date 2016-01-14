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
#ifndef KALMAN_LINEARIZEDMEASUREMENTMODEL_HPP_
#define KALMAN_LINEARIZEDMEASUREMENTMODEL_HPP_

#include "MeasurementModel.hpp"

namespace Kalman {
    template<class StateType>
    class ExtendedKalmanFilter;
    template<class StateType>
    class SquareRootExtendedKalmanFilter;
    
    /**
     * @brief Abstract base class of all linearized (first order taylor expansion) measurement models
     *
     * @param StateType The vector-type of the system state (usually some type derived from Kalman::Vector)
     * @param MeasurementType The vector-type of the measurement (usually some type derived from Kalman::Vector)
     * @param CovarianceBase The class template used for covariance storage (must be either StandardBase or SquareRootBase)
     */
    template<class StateType, class MeasurementType, template<class> class CovarianceBase = StandardBase>
    class LinearizedMeasurementModel : public MeasurementModel<StateType, MeasurementType, CovarianceBase>
    {
        friend class ExtendedKalmanFilter<StateType>;
        friend class SquareRootExtendedKalmanFilter<StateType>;
    public:
        //! Measurement model base
        typedef MeasurementModel<StateType, MeasurementType, CovarianceBase> Base;
        
        //! System state type
        using typename Base::State;
        
        //! Measurement vector type
        using typename Base::Measurement;
        
    protected:
        //! Measurement model jacobian
        Jacobian<Measurement, State> H;
        //! Measurement model noise jacobian
        Jacobian<Measurement, Measurement> V;
        
        /**
         * Callback function for state-dependent update of Jacobi-matrices H and V before each update step
         */
        virtual void updateJacobians( const State& x )
        {
            // No update by default
            (void)x;
        }
    protected:
        LinearizedMeasurementModel()
        {
            H.setIdentity();
            V.setIdentity();
        }
        ~LinearizedMeasurementModel() {}
    };
}

#endif
