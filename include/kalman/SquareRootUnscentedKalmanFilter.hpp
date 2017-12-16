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
#ifndef KALMAN_SQUAREROOTUNSCENTEDKALMANFILTER_HPP_
#define KALMAN_SQUAREROOTUNSCENTEDKALMANFILTER_HPP_

#include "UnscentedKalmanFilterBase.hpp"
#include "SquareRootFilterBase.hpp"

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
    class SquareRootUnscentedKalmanFilter : public UnscentedKalmanFilterBase<StateType>,
                                            public SquareRootFilterBase<StateType>
    {
    public:
        //! Unscented Kalman Filter base type
        typedef UnscentedKalmanFilterBase<StateType> UnscentedBase;
        
        //! Square Root Filter base type
        typedef SquareRootFilterBase<StateType> SquareRootBase;
        
        //! Numeric Scalar Type inherited from base
        using typename UnscentedBase::T;
        
        //! State Type inherited from base
        using typename UnscentedBase::State;
        
        //! Measurement Model Type
        template<class Measurement, template<class> class CovarianceBase>
        using MeasurementModelType = typename UnscentedBase::template MeasurementModelType<Measurement, CovarianceBase>;

        //! System Model Type
        template<class Control, template<class> class CovarianceBase>
        using SystemModelType = typename UnscentedBase::template SystemModelType<Control, CovarianceBase>;
        
    protected:
        //! The number of sigma points (depending on state dimensionality)
        using UnscentedBase::SigmaPointCount;
        
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
        
        //! Square Root of State Covariance
        using SquareRootBase::S;
        
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
        SquareRootUnscentedKalmanFilter(T alpha = T(1), T beta = T(2), T kappa = T(0))
            : UnscentedKalmanFilterBase<StateType>(alpha, beta, kappa)
        {
            // Init covariance to identity
            S.setIdentity();
        }
       
        /**
         * @brief Perform filter prediction step using system model and no control input (i.e. \f$ u = 0 \f$)
         *
         * @param [in] s The System model
         * @return The updated state estimate
         */
        template<class Control, template<class> class CovarianceBase>
        const State& predict( const SystemModelType<Control, CovarianceBase>& s )
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
        const State& predict( const SystemModelType<Control, CovarianceBase>& s, const Control& u )
        {
            // Compute sigma points
            computeSigmaPoints();
            
            // Compute predicted state
            x = this->template computeStatePrediction<Control, CovarianceBase>(s, u);
            
            // Compute predicted covariance
            if(!computeCovarianceSquareRootFromSigmaPoints(x, sigmaStatePoints, s.getCovarianceSquareRoot(), S))
            {
                // TODO: handle numerical error
                assert(false);
            }
            
            // Return predicted state
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
        const State& update( const MeasurementModelType<Measurement, CovarianceBase>& m, const Measurement& z )
        {
            SigmaPoints<Measurement> sigmaMeasurementPoints;
            
            // Predict measurement (and corresponding sigma points)
            Measurement y = this->template computeMeasurementPrediction<Measurement, CovarianceBase>(m, sigmaMeasurementPoints);
            
            // Compute square root innovation covariance
            CovarianceSquareRoot<Measurement> S_y;
            if(!computeCovarianceSquareRootFromSigmaPoints(y, sigmaMeasurementPoints, m.getCovarianceSquareRoot(), S_y))
            {
                // TODO: handle numerical error
                assert(false);
            }
            
            KalmanGain<Measurement> K;
            computeKalmanGain(y, sigmaMeasurementPoints, S_y, K);
            
            // Update state
            x += K * ( z - y );
            
            // Update state covariance
            if(!updateStateCovariance<Measurement>(K, S_y))
            {
                // TODO: handle numerical error
                assert(false);
            }
            
            return this->getState();
        }
        
    protected:
        /**
         * @brief Compute sigma points from current state estimate and state covariance
         * 
         * @note This covers equations (17) and (22) of Algorithm 3.1 in the Paper
         */
        bool computeSigmaPoints()
        {
            // Get square root of covariance
            Matrix<T, State::RowsAtCompileTime, State::RowsAtCompileTime> _S  = S.matrixL().toDenseMatrix();
            
            // Set left "block" (first column)
            sigmaStatePoints.template leftCols<1>() = x;
            // Set center block with x + gamma * S
            sigmaStatePoints.template block<State::RowsAtCompileTime, State::RowsAtCompileTime>(0,1)
                    = ( this->gamma * _S).colwise() + x;
            // Set right block with x - gamma * S
            sigmaStatePoints.template rightCols<State::RowsAtCompileTime>()
                    = (-this->gamma * _S).colwise() + x;
            
            return true;
        }
        
        /**
         * @brief Compute the Covariance Square root from sigma points and noise covariance
         * 
         * @note This covers equations (20) and (21) as well as (25) and (26) of Algorithm 3.1 in the Paper
         * 
         * @param [in] mean The mean predicted state or measurement
         * @param [in] sigmaPoints the predicted sigma state or measurement points
         * @param [in] noiseCov The system or measurement noise covariance (as square root)
         * @param [out] cov The propagated state or innovation covariance (as square root)
         *
         * @return True on success, false if a numerical error is encountered when updating the matrix
         */
        template<class Type>
        bool computeCovarianceSquareRootFromSigmaPoints(const Type& mean, const SigmaPoints<Type>& sigmaPoints, 
                                                        const CovarianceSquareRoot<Type>& noiseCov, CovarianceSquareRoot<Type>& cov)
        {
            // Compute QR decomposition of (transposed) augmented matrix
            Matrix<T, 2*State::RowsAtCompileTime + Type::RowsAtCompileTime, Type::RowsAtCompileTime> tmp;
            tmp.template topRows<2*State::RowsAtCompileTime>() = std::sqrt(this->sigmaWeights_c[1]) * ( sigmaPoints.template rightCols<SigmaPointCount-1>().colwise() - mean).transpose();
            tmp.template bottomRows<Type::RowsAtCompileTime>() = noiseCov.matrixU().toDenseMatrix();

            // TODO: Use ColPivHouseholderQR
            Eigen::HouseholderQR<decltype(tmp)> qr( tmp );
            
            // Set R matrix as upper triangular square root
            cov.setU(qr.matrixQR().template topRightCorner<Type::RowsAtCompileTime, Type::RowsAtCompileTime>());
            
            // Perform additional rank 1 update
            // TODO: According to the paper (Section 3, "Cholesky factor updating") the update
            //       is defined using the square root of the scalar, however the correct result
            //       is obtained when using the weight directly rather than using the square root
            //       It should be checked whether this is a bug in Eigen or in the Paper
            // T nu = std::copysign( std::sqrt(std::abs(sigmaWeights_c[0])), sigmaWeights_c[0]);
            T nu = this->sigmaWeights_c[0];
            cov.rankUpdate( sigmaPoints.template leftCols<1>() - mean, nu );
            
            return (cov.info() == Eigen::Success);
        }
        
        /**
         * @brief Compute the Kalman Gain from predicted measurement and sigma points and the innovation covariance.
         * 
         * @note This covers equations (27) and (28) of Algorithm 3.1 in the Paper
         *
         * We need to solve the equation \f$ K (S_y S_y^T) = P \f$ for \f$ K \f$ using backsubstitution.
         * In order to apply standard backsubstitution using the `solve` method we must bring the
         * equation into the form
         * \f[ AX = B \qquad \text{with } A = S_yS_y^T \f]
         * Thus we transpose the whole equation to obtain
         * \f{align*}{
         *   ( K (S_yS_y^T))^T &= P^T \Leftrightarrow \\
         *   (S_yS_y^T)^T K^T &= P^T \Leftrightarrow \\
         *   (S_yS_y^T) K^T &= P^T
         *\f}
         * Hence our \f$ X := K^T\f$ and \f$ B:= P^T \f$
         *
         * @param [in] y The predicted measurement
         * @param [in] sigmaMeasurementPoints The predicted sigma measurement points
         * @param [in] S_y The innovation covariance as square-root
         * @param [out] K The computed Kalman Gain matrix \f$ K \f$
         */
        template<class Measurement>
        bool computeKalmanGain( const Measurement& y,
                                const SigmaPoints<Measurement>& sigmaMeasurementPoints,
                                const CovarianceSquareRoot<Measurement>& S_y,
                                KalmanGain<Measurement>& K)
        {
            // Note: The intermediate eval() is needed here (for now) due to a bug in Eigen that occurs
            // when Measurement::RowsAtCompileTime == 1 AND State::RowsAtCompileTime >= 8
            decltype(sigmaStatePoints) W = this->sigmaWeights_c.transpose().template replicate<State::RowsAtCompileTime,1>();
            Matrix<T, State::RowsAtCompileTime, Measurement::RowsAtCompileTime> P
                    = (sigmaStatePoints.colwise() - x).cwiseProduct( W ).eval()
                    * (sigmaMeasurementPoints.colwise() - y).transpose();
            
            K = S_y.solve(P.transpose()).transpose();
            return true;
        }
        
        /**
         * @brief Update the state covariance matrix using the Kalman Gain and the Innovation Covariance
         * 
         * @note This covers equations (29) and (30) of Algorithm 3.1 in the Paper
         *
         * @param [in] K The computed Kalman Gain matrix
         * @param [in] S_y The innovation covariance as square-root
         * @return True on success, false if a numerical error is encountered when updating the matrix
         */
        template<class Measurement>
        bool updateStateCovariance(const KalmanGain<Measurement>& K, const CovarianceSquareRoot<Measurement>& S_y)
        {
            KalmanGain<Measurement> U = K * S_y.matrixL();
            for(int i = 0; i < U.cols(); ++i)
            {
                S.rankUpdate( U.col(i), -1 );
                
                if( S.info() == Eigen::NumericalIssue )
                {
                    // TODO: handle numerical issues in some sensible way
                    return false;
                }
            }
            
            return true;
        }
    };
}

#endif
