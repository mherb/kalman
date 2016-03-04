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
#ifndef KALMAN_SQUAREROOTEXTENDEDKALMANFILTER_HPP_
#define KALMAN_SQUAREROOTEXTENDEDKALMANFILTER_HPP_

#include "KalmanFilterBase.hpp"
#include "SquareRootFilterBase.hpp"
#include "LinearizedSystemModel.hpp"
#include "LinearizedMeasurementModel.hpp"

namespace Kalman {
    
    /**
     * @brief Square-Root Extended Kalman Filter (SR-EKF)
     * 
     * This implementation is based upon [An Introduction to the Kalman Filter](https://www.cs.unc.edu/~welch/media/pdf/kalman_intro.pdf)
     * by Greg Welch and Gary Bishop.
     *
     * @param StateType The vector-type of the system state (usually some type derived from Kalman::Vector)
     */
    template<class StateType>
    class SquareRootExtendedKalmanFilter : public KalmanFilterBase<StateType>,
                                           public SquareRootFilterBase<StateType>
    {
    public:
        //! Kalman Filter base type
        typedef KalmanFilterBase<StateType> KalmanBase;
        //! SquareRoot Filter base type
        typedef SquareRootFilterBase<StateType> SquareRootBase;
        
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
        using KalmanGain = Kalman::KalmanGain<State, Measurement>;
        
    protected:
        //! State estimate
        using KalmanBase::x;
        //! State covariance square root
        using SquareRootBase::S;
        
    public:
        /**
         * @brief Constructor
         */
        SquareRootExtendedKalmanFilter()
        {
            // Setup covariance
            S.setIdentity();
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
            computePredictedCovarianceSquareRoot<State>(s.F, S, s.W, s.getCovarianceSquareRoot(), S);
            
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
            CovarianceSquareRoot<Measurement> S_y;
            computePredictedCovarianceSquareRoot<Measurement>(m.H, S, m.V, m.getCovarianceSquareRoot(), S_y);
            
            // compute kalman gain
            KalmanGain<Measurement> K;
            computeKalmanGain<Measurement>(m.H, S_y, K);
            
            // UPDATE STATE ESTIMATE AND COVARIANCE
            // Update state using computed kalman gain and innovation
            x += K * ( z - m.h( x ) );
            
            // Update covariance
            updateStateCovariance<Measurement>(K, m.H);
            
            // return updated state estimate
            return this->getState();
        }
    protected:
        
        /**
         * @brief Compute the predicted state or innovation covariance (as square root)
         *
         * The formula for computing the propagated square root covariance matrix can be
         * deduced in a very straight forward way using the method outlined in the
         * [Unscented Kalman Filter Tutorial](https://cse.sc.edu/~terejanu/files/tutorialUKF.pdf)
         * by Gabriel A. Terejanu for the UKF (Section 3, formulas (27) and (28)).
         * 
         * Starting from the standard update formula
         * 
         *     \f[ \hat{P} = FPF^T + WQW^T \f]
         * 
         * and using the square-root decomposition \f$ P = SS^T \f$ with \f$S\f$ being lower-triagonal
         * as well as the (lower triangular) square root \f$\sqrt{Q}\f$ of \f$Q\f$ this can be formulated as
         *
         *     \f{align*}{
         *         \hat{P}  &= FSS^TF^T + W\sqrt{Q}\sqrt{Q}^TW^T \\
         *                  &= FS(FS)^T + W\sqrt{Q}(W\sqrt{Q})^T \\
         *                  &=  \begin{bmatrix} FS & W\sqrt{Q} \end{bmatrix}
         *                      \begin{bmatrix} (FS)^T \\ (W\sqrt{Q})^T \end{bmatrix}
         *     \f}
         *
         * The blockmatrix
         *
         * \f[ \begin{bmatrix} (FS)^T \\ (W\sqrt{Q})^T \end{bmatrix} \in \mathbb{R}^{2n \times n} \f]
         *
         * can then be decomposed into a product of matrices \f$OR\f$ with \f$O\f$ being orthogonal
         * and \f$R\f$ being upper triangular (also known as QR decompositon). Using this \f$\hat{P}\f$
         * can be written as
         * 
         *     \f{align*}{
         *         \hat{P}  &=  \begin{bmatrix} FS & W\sqrt{Q} \end{bmatrix}
         *                      \begin{bmatrix} (FS)^T \\ (W\sqrt{Q})^T \end{bmatrix} \\
         *                  &= (OR)^TOR \\
         *                  &= R^T \underbrace{O^T O}_{=I}R \\
         *                  &= LL^T \qquad \text{ with } L := R^T
         *     \f}
         *
         * Thus the lower-triangular square root of \f$\hat{P}\f$ is equaivalent to the transpose
         * of the upper-triangular matrix obtained from QR-decompositon of the augmented block matrix
         *
         *     \f[ \begin{bmatrix} FS & W\sqrt{Q} \end{bmatrix}^T = \begin{bmatrix} S^T F^T \\ \sqrt{Q}^T W^T \end{bmatrix} \in \mathbb{R}^{2n \times n} \f]
         * 
         * The same can be applied for the innovation covariance by replacing the jacobians and
         * the noise covariance accordingly.
         * 
         * @param [in] A The jacobian of state transition or measurement function w.r.t. state or measurement
         * @param [in] S The state covariance (as square root)
         * @param [in] B The jacobian of state transition or measurement function w.r.t. state or measurement
         * @param [in] R The system model or measurement noise covariance (as square root)
         * @param [out] S_pred The predicted covariance (as square root)
         * @return True on success, false on failure due to numerical issue
         */
        template<class Type>
        bool computePredictedCovarianceSquareRoot(  const Jacobian<Type, State>& A, const CovarianceSquareRoot<State>& S,
                                                    const Jacobian<Type, Type>& B,  const CovarianceSquareRoot<Type>& R,
                                                    CovarianceSquareRoot<Type>& S_pred)
        {
            // Compute QR decomposition of (transposed) augmented matrix
            Matrix<T, State::RowsAtCompileTime + Type::RowsAtCompileTime, Type::RowsAtCompileTime> tmp;
            tmp.template topRows<State::RowsAtCompileTime>()   = S.matrixU() * A.transpose();
            tmp.template bottomRows<Type::RowsAtCompileTime>() = R.matrixU() * B.transpose();
            
            // TODO: Use ColPivHouseholderQR
            Eigen::HouseholderQR<decltype(tmp)> qr( tmp );
            
            // Set S_pred matrix as upper triangular square root
            S_pred.setU(qr.matrixQR().template topRightCorner<Type::RowsAtCompileTime, Type::RowsAtCompileTime>());
            return true;
        }
        
        /**
         * @brief Compute Kalman gain from state covariance, innovation covariance and measurement jacobian
         * 
         * The formula for the kalman gain using square-root covariances can be deduced from the
         * original update formula
         *
         *      \f[ K = P H^T P_{yy}^{-1} \Leftrightarrow K P_{yy} = P H^T \f]
         * 
         * with \f$P_yy\f$ being the innovation covariance. With \f$P = SS^T\f$ and \f$P_yy = S_yS_y^T\f$
         * and transposition of the whole equation the above can be formulated as
         *
         *  \f{align*}{
         *      K (S_yS_y^T) &= SS^TH^T \\
         *      (K(S_yS_y^T))^T &= (SS^TH^T)^T \\
         *      (S_yS_y^T)^TK^T &= H(SS^T)^T \\
         *      \underbrace{S_yS_y^T}_{=P_{yy}} K^T &= HSS^T
         *  \f}
         * 
         * Since \f$S_y\f$ is lower triangular, the above can be solved using backsubstitution
         * in an efficient manner.
         *
         * @param [in] H The jacobian of the measurement function w.r.t. the state
         * @param [in] S_y The Innovation covariance square root
         * @param [out] K The computed Kalman Gain
         */
        template<class Measurement>
        bool computeKalmanGain( const Jacobian<Measurement, State>& H,
                                const CovarianceSquareRoot<Measurement>& S_y,
                                KalmanGain<Measurement>& K)
        {
            // Solve using backsubstitution
            // AX=B with B = HSS^T and X = K^T and A = S_yS_y^T
            K = S_y.solve(H*S.reconstructedMatrix()).transpose();
            return true;
        }
        
        /**
         * @brief Update state covariance using Kalman gain
         *
         * @param [in] K The Kalman gain
         * @param [in] H The jacobian of the measurement function w.r.t. the state
         */
        template<class Measurement>
        bool updateStateCovariance( const KalmanGain<Measurement>& K,
                                    const Jacobian<Measurement, State>& H)
        {
            // TODO: update covariance without using decomposition
            Matrix<T, State::RowsAtCompileTime, State::RowsAtCompileTime> P = S.reconstructedMatrix();
            S.compute( (P - K * H * P).eval() );
            return (S.info() == Eigen::Success);
        }
    };
}

#endif
