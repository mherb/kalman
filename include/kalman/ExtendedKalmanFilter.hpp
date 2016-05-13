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
#ifndef KALMAN_EXTENDEDKALMANFILTER_HPP_
#define KALMAN_EXTENDEDKALMANFILTER_HPP_

#include "KalmanFilterBase.hpp"
#include "StandardFilterBase.hpp"
#include "LinearizedSystemModel.hpp"
#include "LinearizedMeasurementModel.hpp"
#include "Model.hpp"
#include "AutoDiff.hpp"

namespace Kalman {
    
    /**
     * @brief Extended Kalman Filter (EKF)
     * 
     * This implementation is based upon [An Introduction to the Kalman Filter](https://www.cs.unc.edu/~welch/media/pdf/kalman_intro.pdf)
     * by Greg Welch and Gary Bishop.
     *
     * @param StateType The vector-type of the system state (usually some type derived from Kalman::Vector)
     */
    template<class StateType>
    class ExtendedKalmanFilter : public KalmanFilterBase<StateType>,
                                 public StandardFilterBase<StateType>
    {
    public:
        //! Kalman Filter base type
        typedef KalmanFilterBase<StateType> KalmanBase;
        //! Standard Filter base type
        typedef StandardFilterBase<StateType> StandardBase;
        
        //! Numeric Scalar Type inherited from base
        using typename KalmanBase::T;
        
        //! State Type inherited from base
        using typename KalmanBase::State;
        
        //! Linearized Measurement Model Type
        template<class Measurement, template<class> class CovarianceBase>
        using MeasurementModelType = LinearizedMeasurementModel<State, Measurement, CovarianceBase>;
        
        //! Linearized System Model Type
        template<class Control, template<class> class CovarianceBase>
        using SystemModelType = LinearizedSystemModel<State, Control, CovarianceBase>;
        
    protected:
        //! Kalman Gain Matrix Type
        template<class Measurement>
        using KalmanGain = Kalman::KalmanGain<T, State, Measurement>;

        //! Covariance type alias
        template<typename Type>
        using Covariance = Kalman::Covariance<T, Type>;

        //! AutoDiff Helper for prediction step
        using PredictAutoDiff = AutoDiffJacobian<T, State>;

        //! AutoDiff Helper for update step
        template<class Measurement>
        using UpdateAutoDiff = AutoDiffJacobian<T, State, Measurement>;

        //! State Jacobian type alias
        typedef Jacobian<T, State> StateJacobian;

        //! Measurement Jacobian type alias
        template<class Measurement>
        using MeasurementJacobian = Jacobian<T, State, Measurement>;

    protected:
        //! State Estimate
        using KalmanBase::x;
        //! State Covariance Matrix
        using StandardBase::P;
        
    public:
        /**
         * @brief Constructor
         */
        ExtendedKalmanFilter()
        {
            // Setup state and covariance
            P.setIdentity();
        }
        
        /**
         * @brief Perform filter prediction step using system model and no control input (i.e. \f$ u = 0 \f$)
         *
         * @param [in] s The System model
         * @return The updated state estimate
         */
        template<class Control, template<class> class CovarianceBase>
        const State& predict( SystemModelType<Control, CovarianceBase>& s )
        {
            // predict state (without control)
            Control u;
            u.setZero();
            return predict( s, u );
        }
        
        /**
         * @brief Perform filter prediction step using control input \f$u\f$ and corresponding system model
         *
         * @param [in] s The System model
         * @param [in] u The Control input vector
         * @return The updated state estimate
         */
        template<class Control, template<class> class CovarianceBase>
        const State& predict( SystemModelType<Control, CovarianceBase>& s, const Control& u )
        {
            s.updateJacobians( x, u );
            
            // predict state
            x = s.f(x, u);
            
            // predict covariance
            P  = ( s.F * P * s.F.transpose() ) + ( s.W * s.getCovariance() * s.W.transpose() );
            
            // return state prediction
            return this->getState();
        }
        
        /**
         * @brief Perform filter update step using measurement \f$z\f$ and corresponding measurement model
         *
         * @param [in] m The Measurement model
         * @param [in] z The measurement vector
         * @return The updated state estimate
         */
        template<class Measurement, template<class> class CovarianceBase>
        const State& update( MeasurementModelType<Measurement, CovarianceBase>& m, const Measurement& z )
        {
            m.updateJacobians( x );
            
            // COMPUTE KALMAN GAIN
            // compute innovation covariance
            Covariance<Measurement> S = ( m.H * P * m.H.transpose() ) + ( m.V * m.getCovariance() * m.V.transpose() );
            
            // compute kalman gain
            KalmanGain<Measurement> K = P * m.H.transpose() * S.inverse();
            
            // UPDATE STATE ESTIMATE AND COVARIANCE
            // Update state using computed kalman gain and innovation
            x += K * ( z - m.h( x ) );
            
            // Update covariance
            P -= K * m.H * P;
            
            // return updated state estimate
            return this->getState();
        }

        /**
         * @brief Perform filter prediction step using system model and additional inputs
         *
         * @param [in] system The System model
         * @param [in] args Additional model prediction arguments (i.e. control input, time-step, ...)
         */
        template<class System, typename... Args>
        typename std::enable_if< Model::isTemplateModel<System, State, Args...>::value >::type
        predict(System& system, Args&&... args)
        {
            State prediction;
            StateJacobian F;

            // Compute prediction and jacobian
            computePrediction(prediction, F, system, std::forward<Args>(args)...);

            // Update state
            x = prediction;

            // Predict Covariance
            // Note: For simplicity, this assumes purely additive noise at the moment.
            // For non-additive noise, one needs another Jacobian W of system.predict
            // w.r.t. the noise vector w. This may be added in future versions.
            P  = ( F * P * F.transpose() ) + Model::template getCovariance<System, T, State>(system);
        }

        /**
         * @brief Perform filter update step using measurement \f$z\f$ and corresponding measurement model
         *
         * @param [in] model The Measurement model
         * @param [in] z The measurement vector
         * @param [in] args Additional measurement arguments
         */
        template<class MeasurementModel, class Measurement, typename... Args>
        typename std::enable_if< Model::isTemplateModel<MeasurementModel, State, Measurement, Args...>::value >::type
        update(MeasurementModel& model, const Measurement& z, Args&&... args)
        {
            // Compute measurement prediction and jacobian
            Measurement prediction;
            MeasurementJacobian<Measurement> H;
            computeMeasurement(prediction, H, model, std::forward<Args>(args)...);

            // COMPUTE KALMAN GAIN
            // compute innovation covariance
            Covariance<Measurement> S = ( H * P * H.transpose() )
                                        + Model::template getCovariance<MeasurementModel, T, Measurement>(model);

            // compute kalman gain
            KalmanGain<Measurement> K = P * H.transpose() * S.inverse();

            // UPDATE STATE ESTIMATE AND COVARIANCE
            // Update state using computed kalman gain and innovation
            x += K * ( z - prediction ).template cast<T>();

            // Update covariance
            P -= K * H * P;
        }
    protected:
        /**
         * @brief Compute Prediction using Auto-Diff Jacobian
         * @param [out] prediction The predicted state
         * @param [out] jacobian Jacobian of prediction function w.r.t. state
         * @param [in] system The system model instance
         * @param [in] args Additional prediction arguments (i.e. control input)
         */
        template<class System, typename... Args>
        typename std::enable_if< !Model::hasJacobian<System, State, Args...>::value >::type
        computePrediction(State& prediction, StateJacobian& jacobian, System& system, Args&&... args)
        {
            typename PredictAutoDiff::Function func
                    = [&](const typename PredictAutoDiff::ActiveInput& in,
                          typename PredictAutoDiff::ActiveOutput& out)
                    {
                        system.predict(in, std::forward<Args>(args)..., out);
                    };

            // Note: This doesn't work since the template types cannot be inferred
            // using namespace std::placeholders;
            // typename PredictAutoDiff::Function func
            //         = std::bind(&System::predict, &system, _1, std::ref(args)..., _2);

            // Compute prediction (with jacobian)
            PredictAutoDiff predictAutoDiff(func);
            predictAutoDiff(x, prediction, jacobian);
        }

        /**
         * @brief Compute Prediction using explicit Jacobian
         * @param [out] prediction The predicted state
         * @param [out] jacobian Jacobian of prediction function w.r.t. state
         * @param [in] system The system model instance
         * @param [in] args Additional prediction arguments (i.e. control input)
         */
        template<class System, typename... Args>
        typename std::enable_if< Model::hasJacobian<System, State, Args...>::value >::type
        computePrediction(State& prediction, StateJacobian& jacobian, System& system, Args&&... args)
        {
            system.predict(x, std::forward<Args>(args)..., prediction);
            jacobian = system.getJacobian(x, std::forward<Args>(args)...);
        }

        /**
         * @brief Compute Measurement prediction using Auto-Diff Jacobian
         * @param [out] prediction The predicted/expected measurement vector
         * @param [out] jacobian Jacobian of measurement function w.r.t state
         * @param [in] model Measurement model instance
         * @param [in] args Additional measurement function arguments
         */
        template<class MeasurementModel, class Measurement, typename... Args>
        typename std::enable_if< !Model::hasJacobian<MeasurementModel, Measurement, Args...>::value >::type
        computeMeasurement(Measurement& prediction, MeasurementJacobian<Measurement>& jacobian, MeasurementModel& model, Args&&... args)
        {
            typename UpdateAutoDiff<Measurement>::Function func
                    = [&](const typename UpdateAutoDiff<Measurement>::ActiveInput& in,
                          typename UpdateAutoDiff<Measurement>::ActiveOutput& out)
                    {
                        model.measure(in, std::forward(args)..., out);
                    };

            // Compute measurement-prediction (with jacobian)
            UpdateAutoDiff<Measurement> updateAutoDiff(func);
            updateAutoDiff(x, prediction, jacobian);
        }

        /**
         * @brief Compute Measurement prediction using explicit Jacobian
         * @param [out] prediction The predicted/expected measurement vector
         * @param [out] jacobian Jacobian of measurement function w.r.t state
         * @param [in] model Measurement model instance
         * @param [in] args Additional measurement function arguments
         */
        template<class MeasurementModel, class Measurement, typename... Args>
        typename std::enable_if< Model::hasJacobian<MeasurementModel, Measurement, Args...>::value >::type
        computeMeasurement(Measurement& prediction, MeasurementJacobian<Measurement>& jacobian, MeasurementModel& model, Args&&... args)
        {
            model.measure(x, std::forward(args)..., prediction);
            jacobian = model.getJacobian(x, std::forward(args)...);
        }
    };
}

#endif
