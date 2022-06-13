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
#ifndef KALMAN_MAHALANOBIS_UNSCENTED_KALMAN_FILTER_HPP_
#define KALMAN_MAHALANOBIS_UNSCENTED_KALMAN_FILTER_HPP_

#include "UnscentedKalmanFilter.hpp"
#include "MahalanobisMeasurementModel.hpp"
#include <memory>

namespace Kalman {

    /**
     * @brief Unscented Kalman Filter (UKF)
     *
     * @note It is recommended to use the square-root implementation SquareRootUnscentedKalmanFilter of this filter
     *
     * This implementation is based upon [The square-root unscented Kalman filter for state and parameter-estimation](http://dx.doi.org/10.1109/ICASSP.2001.940586) by Rudolph van der Merwe and Eric A. Wan.
     * Whenever "the paper" is referenced within this file then this paper is meant.
     *
     * @param StateType The vector-type of the system state (usually some type derived from Kalman::Vector)
     */
    template<class StateType>
    class UnscentedMahalanobisKalmanFilter : public UnscentedKalmanFilter<StateType>
    {
    public:
        //! Unscented Kalman Filter type
        typedef UnscentedKalmanFilter<StateType> UnscentedFilterType;

        //! Unscented Kalman Filter base type
        typedef UnscentedKalmanFilterBase<StateType> UnscentedBase;

        //! Numeric Scalar Type inherited from base
        using typename UnscentedBase::T;

        //! State Type inherited from base
        using typename UnscentedBase::State;

        //! Robust Measurement Model Type
        template<class Measurement, template<class> class CovarianceBase>
        using MahalanobisMeasurementModelType = MahalanobisMeasurementModel<State, Measurement, CovarianceBase>;

        //! Measurement Model Type
        template<class Measurement, template<class> class CovarianceBase>
        using MeasurementModelType = MeasurementModel<State, Measurement, CovarianceBase>;

    protected:
        //! Matrix type containing the sigma state or measurement points
        template<class Type>
        using SigmaPoints = typename UnscentedBase::template SigmaPoints<Type>;

        //! Kalman Gain Matrix Type
        template<class Measurement>
        using KalmanGain = Kalman::KalmanGain<State, Measurement>;

    protected:
        //! State Estimate (member)
        using UnscentedBase::x;
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
        UnscentedMahalanobisKalmanFilter(T alpha = T(1), T beta = T(2), T kappa = T(0)):
        UnscentedFilterType(alpha, beta, kappa)
        {}

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
            return UnscentedFilterType::template update(m, z);
        }

        /**
         * @brief Perform filter update (with outlier prevention) step using measurement \f$z\f$ and corresponding measurement model
         *
         * @param [in] m    The Measurement model
         * @param [in] z    The measurement vector
         * @param [out] accepted    flag indicating wether or not the measurement was deemed to be an outlier,
         *                          and hence _not_ used for filter update
         * @return The (possibly, depending on the accepted flag) updated state estimate
         */
        template<class Measurement, template<class> class CovarianceBase>
        const State& update(
            const MahalanobisMeasurementModelType<Measurement, CovarianceBase>& m,
            const Measurement& z,
            bool& accepted
        )
        {
            SigmaPoints<Measurement> sigmaMeasurementPoints;
            const Measurement y{
                UnscentedFilterType::template computeMeasurementPrediction<Measurement, CovarianceBase>(m, sigmaMeasurementPoints)
            };

            const Measurement innovation{z - y};

            Covariance<Measurement> innovation_covariance;
            UnscentedFilterType::template computeCovarianceFromSigmaPoints(y, sigmaMeasurementPoints, m.getCovariance(), innovation_covariance);
            const T squared_mahalanobis_distance{innovation.transpose() * innovation_covariance.inverse() * innovation};

            accepted = squared_mahalanobis_distance <= m.getChiSquaredThreshold();
            if (accepted)
            {
                KalmanGain<Measurement> K;
                UnscentedFilterType::template computeKalmanGain(innovation, sigmaMeasurementPoints, innovation_covariance, K);

                x += K * innovation;
                UnscentedFilterType::template updateStateCovariance<Measurement>(K, innovation_covariance);
            }
            return x;
        }
    };
}

#endif
