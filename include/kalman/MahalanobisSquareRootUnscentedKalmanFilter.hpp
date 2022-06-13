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
#ifndef KALMAN_MAHALANOBIS_SQUARE_ROOT_UNSCENTED_KALMAN_FILTER_HPP_
#define KALMAN_MAHALANOBIS_SQUARE_ROOT_UNSCENTED_KALMAN_FILTER_HPP_

#include "SquareRootUnscentedKalmanFilter.hpp"
#include "MahalanobisMeasurementModel.hpp"
#include <iostream>
namespace Kalman {

    /**
     * @brief Square Root Unscented Kalman Filter (SR-UKF)
     *
     * @note This is the square-root implementation variant of UnscentedKalmanFilter
     *
     * This implementation is based upon [The square-root unscented Kalman filter for state and parameter-estimation](http://dx.doi.org/10.1109/ICASSP.2001.940586) by Rudolph van der Merwe and Eric A. Wan.
     * Whenever "the paper" is referenced within this file then this paper is meant.
     *
     * @param StateType The vector-type of the system state (usually some type derived from Kalman::Vector)
     */
    template<class StateType>
    class MahalanobisSquareRootUnscentedKalmanFilter : public SquareRootUnscentedKalmanFilter<StateType>
    {
    public:
        //! Unscented Kalman Filter type
        typedef SquareRootUnscentedKalmanFilter<StateType> SquareRootUnscentedKalmanFilterType;

        //! Unscented Kalman Filter base type
        typedef UnscentedKalmanFilterBase<StateType> UnscentedBase;

        //! Square Root Filter base type
        typedef SquareRootFilterBase<StateType> SquareRootBase;

        //! Numeric Scalar Type inherited from base
        using typename UnscentedBase::T;

        //! State Type inherited from base
        using typename UnscentedBase::State;

        //! Robust Measurement Model Type
        template<class Measurement, template<class> class CovarianceBase>
        using MahalanobisMeasurementModelType = MahalanobisMeasurementModel<State, Measurement, CovarianceBase>;

        //! Measurement Model Type
        template<class Measurement, template<class> class CovarianceBase>
        using MeasurementModelType = typename UnscentedBase::template MeasurementModelType<Measurement, CovarianceBase>;

    protected:
        //! Matrix type containing the sigma state or measurement points
        template<class Type>
        using SigmaPoints = typename UnscentedBase::template SigmaPoints<Type>;

        //! Kalman Gain Matrix Type
        template<class Measurement>
        using KalmanGain = Kalman::KalmanGain<State, Measurement>;

    protected:
        // Member variables

        //! State Estimate
        using UnscentedBase::x;

        //! Sigma points (state)
        using UnscentedBase::sigmaStatePoints;

    public:
        /**
         * Constructor
         *
         * See paper for detailed parameter explanation
         *
         * @param alpha Scaling parameter for spread of sigma points (usually \f$ 1E-4 \leq \alpha \leq 1 \f$)
         * @param beta Parameter for prior knowledge about the distribution (\f$ \beta = 2 \f$ is optimal for Gaussian)
         * @param kappa Secondary scaling parameter (usually 0)
         */
        MahalanobisSquareRootUnscentedKalmanFilter(T alpha = T(1), T beta = T(2), T kappa = T(0)):
        SquareRootUnscentedKalmanFilter<StateType>{alpha, beta, kappa} {}

        /**
         * @brief Perform filter update step using measurement \f$z\f$ and corresponding measurement model
         *
         * @param [in] m The Measurement model
         * @param [in] z The measurement vector
         * @return The updated state estimate
         */
        template<class Measurement, template<class> class CovarianceBase>
        const State& update( const MeasurementModelType<Measurement, CovarianceBase>& m, const Measurement& z )
        {
            return SquareRootUnscentedKalmanFilterType::template update(m, z);
        }

                /**
         * @brief Perform filter update step using measurement \f$z\f$ and corresponding measurement model
         *
         * @param [in] m    The Measurement model
         * @param [in] z    The measurement vector
         * @param [out] accepted    flag indicating wether or not the measurement was deemed to be an outlier,
         *                          and hence _not_ used for filter update
         * @return The (possibly, depending on the accepted flag) updated state estimate
         */
        template<class Measurement, template<class> class CovarianceBase>
        const State& update(const MahalanobisMeasurementModelType<Measurement, CovarianceBase>& m, const Measurement& z, bool& accepted)
        {
            SigmaPoints<Measurement> sigmaMeasurementPoints;
            const Measurement y{
                SquareRootUnscentedKalmanFilterType::template computeMeasurementPrediction<Measurement, CovarianceBase>(m, sigmaMeasurementPoints)
            };

            const Measurement innovation{z - y};

            CovarianceSquareRoot<Measurement> innovation_covariance_sqrt;
            SquareRootUnscentedKalmanFilterType::template computeCovarianceSquareRootFromSigmaPoints(y, sigmaMeasurementPoints, m.getCovarianceSquareRoot(), innovation_covariance_sqrt);

            const T squared_mahalanobis_distance{innovation.transpose() * innovation_covariance_sqrt.reconstructedMatrix().inverse() * innovation};

            accepted = squared_mahalanobis_distance <= m.getChiSquaredThreshold();
            if (accepted)
            {
                KalmanGain<Measurement> K;
                SquareRootUnscentedKalmanFilterType::template computeKalmanGain(innovation, sigmaMeasurementPoints, innovation_covariance_sqrt, K);

                x += K * innovation;
                SquareRootUnscentedKalmanFilterType::template updateStateCovariance<Measurement>(K, innovation_covariance_sqrt);
            }
            return x;
        }
    };
}

#endif
