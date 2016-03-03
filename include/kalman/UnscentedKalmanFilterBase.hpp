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
#ifndef KALMAN_UNSCENTEDKALMANFILTERBASE_HPP_
#define KALMAN_UNSCENTEDKALMANFILTERBASE_HPP_

#include <cassert>

#include "KalmanFilterBase.hpp"
#include "SystemModel.hpp"
#include "MeasurementModel.hpp"

namespace Kalman {
    
    /**
     * @brief Abstract Base for Unscented Kalman Filters
     *
     * This implementation is based upon [The square-root unscented Kalman filter for state and parameter-estimation](http://dx.doi.org/10.1109/ICASSP.2001.940586) by Rudolph van der Merwe and Eric A. Wan.
     * Whenever "the paper" is referenced within this file then this paper is meant.
     * 
     * @param StateType The vector-type of the system state (usually some type derived from Kalman::Vector)
     */
    template<class StateType>
    class UnscentedKalmanFilterBase : public KalmanFilterBase<StateType>
    {
    public:
        // Public typedefs
        //! Kalman Filter base type
        typedef KalmanFilterBase<StateType> Base;
        
        //! Numeric Scalar Type inherited from base
        using typename Base::T;
        
        //! State Type inherited from base
        using typename Base::State;
    
        //! Measurement Model Type
        template<class Measurement, template<class> class CovarianceBase>
        using MeasurementModelType = MeasurementModel<State, Measurement, CovarianceBase>;
        
        //! System Model Type
        template<class Control, template<class> class CovarianceBase>
        using SystemModelType = SystemModel<State, Control, CovarianceBase>;
        
    protected:
        // Protected typedefs
        //! The number of sigma points (depending on state dimensionality)
        static constexpr int SigmaPointCount = 2 * State::RowsAtCompileTime + 1;
        
        //! Vector containg the sigma scaling weights
        typedef Vector<T, SigmaPointCount> SigmaWeights;
        
        //! Matrix type containing the sigma state or measurement points
        template<class Type>
        using SigmaPoints = Matrix<T, Type::RowsAtCompileTime, SigmaPointCount>;

    protected:
        // Member variables
        
        //! State Estimate
        using Base::x;
        
        //! Sigma weights (m)
        SigmaWeights sigmaWeights_m;
        //! Sigma weights (c)
        SigmaWeights sigmaWeights_c;
        
        //! Sigma points (state)
        SigmaPoints<State> sigmaStatePoints;
        
        // Weight parameters
        T alpha;    //!< Scaling parameter for spread of sigma points (usually \f$ 1E-4 \leq \alpha \leq 1 \f$)
        T beta;     //!< Parameter for prior knowledge about the distribution (\f$ \beta = 2 \f$ is optimal for Gaussian)
        T kappa;    //!< Secondary scaling parameter (usually 0)
        T gamma;    //!< \f$ \gamma = \sqrt{L + \lambda} \f$ with \f$ L \f$ being the state dimensionality
        T lambda;   //!< \f$ \lambda = \alpha^2 ( L + \kappa ) - L\f$ with \f$ L \f$ being the state dimensionality
        
    protected:
        /**
         * Constructor
         * 
         * See paper for parameter explanation
         * 
         * @param _alpha Scaling parameter for spread of sigma points (usually 1e-4 <= alpha <= 1)
         * @param _beta Parameter for prior knowledge about the distribution (beta = 2 is optimal for Gaussian)
         * @param _kappa Secondary scaling parameter (usually 0)
         */
        UnscentedKalmanFilterBase(T _alpha = T(1), T _beta = T(2), T _kappa = T(0)) : alpha(_alpha), beta(_beta), kappa(_kappa)
        {
            // Pre-compute all weights
            computeWeights();
            
            // Setup state and covariance
            x.setZero();
        }

        /**
         * @brief Compute predicted state using system model and control input
         *
         * @param [in] s The System Model
         * @param [in] u The Control input
         * @return The predicted state
         */
        template<class Control, template<class> class CovarianceBase>
        State computeStatePrediction(const SystemModelType<Control, CovarianceBase>& s, const Control& u)
        {
            // Pass each sigma point through non-linear state transition function
            computeSigmaPointTransition(s, u);
            
            // Compute predicted state from predicted sigma points
            return computePredictionFromSigmaPoints<State>(sigmaStatePoints);
        }
        
        /**
         * @brief Compute predicted measurement using measurement model and predicted sigma measurements
         *
         * @param [in] m The Measurement Model
         * @param [in] sigmaMeasurementPoints The predicted sigma measurement points
         * @return The predicted measurement
         */
        template<class Measurement, template<class> class CovarianceBase>
        Measurement computeMeasurementPrediction(const MeasurementModelType<Measurement, CovarianceBase>& m, SigmaPoints<Measurement>& sigmaMeasurementPoints)
        {
            // Predict measurements for each sigma point
            computeSigmaPointMeasurements<Measurement>(m, sigmaMeasurementPoints);
            
            // Predict measurement from sigma measurement points
            return computePredictionFromSigmaPoints<Measurement>(sigmaMeasurementPoints);
        }
        
        /**
         * @brief Compute sigma weights
         */
        void computeWeights()
        {
            T L = T(State::RowsAtCompileTime);
            lambda = alpha * alpha * ( L + kappa ) - L;
            gamma = std::sqrt( L + lambda );
            
            // Make sure L != -lambda to avoid division by zero
            assert( std::abs(L + lambda) > 1e-6 );
            
            // Make sure L != -kappa to avoid division by zero
            assert( std::abs(L + kappa) > 1e-6 );
            
            T W_m_0 = lambda / ( L + lambda );
            T W_c_0 = W_m_0 + (T(1) - alpha*alpha + beta);
            T W_i   = T(1) / ( T(2) * alpha*alpha * (L + kappa) );
            
            // Make sure W_i > 0 to avoid square-root of negative number
            assert( W_i > T(0) );
            
            sigmaWeights_m[0] = W_m_0;
            sigmaWeights_c[0] = W_c_0;
            
            for(int i = 1; i < SigmaPointCount; ++i)
            {
                sigmaWeights_m[i] = W_i;
                sigmaWeights_c[i] = W_i;
            }
        }
        
        /**
         * @brief Predict expected sigma states from current sigma states using system model and control input
         * 
         * @note This covers equation (18) of Algorithm 3.1 in the Paper
         *
         * @param [in] s The System Model
         * @param [in] u The Control input
         */
        template<class Control, template<class> class CovarianceBase>
        void computeSigmaPointTransition(const SystemModelType<Control, CovarianceBase>& s, const Control& u)
        {
            for( int i = 0; i < SigmaPointCount; ++i )
            {
                sigmaStatePoints.col(i) = s.f( sigmaStatePoints.col(i), u );
            }
        }
        
        /**
         * @brief Predict the expected sigma measurements from predicted sigma states using measurement model
         * 
         * @note This covers equation (23) of Algorithm 3.1 in the Paper
         *
         * @param [in] m The Measurement model
         * @param [out] sigmaMeasurementPoints The struct of expected sigma measurements to be computed
         */
        template<class Measurement, template<class> class CovarianceBase>
        void computeSigmaPointMeasurements(const MeasurementModelType<Measurement, CovarianceBase>& m, SigmaPoints<Measurement>& sigmaMeasurementPoints)
        {
            for( int i = 0; i < SigmaPointCount; ++i )
            {
                sigmaMeasurementPoints.col(i) = m.h( sigmaStatePoints.col(i) );
            }
        }
        
        /**
         * @brief Compute state or measurement prediciton from sigma points using pre-computed sigma weights
         * 
         * @note This covers equations (19) and (24) of Algorithm 3.1 in the Paper
         *
         * @param [in] sigmaPoints The computed sigma points of the desired type (state or measurement)
         * @return The prediction
         */
        template<class Type>
        Type computePredictionFromSigmaPoints(const SigmaPoints<Type>& sigmaPoints)
        {
            // Use efficient matrix x vector computation
            return sigmaPoints * sigmaWeights_m;
        }
    };
}

#endif
