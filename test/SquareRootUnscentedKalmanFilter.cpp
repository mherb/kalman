#include "TestHelper.h"

#define private public
#define protected public

#include <kalman/SquareRootUnscentedKalmanFilter.hpp>

using namespace Kalman;

typedef float T;

TEST(SquareRootUnscentedKalmanFilter, init) {
    auto ukf = SquareRootUnscentedKalmanFilter<Vector<T, 3>>(1,2,1);
    ASSERT_TRUE(ukf.S.isIdentity()); // S should be identity

    // Same as above, but with general matrix type instead of vector
    auto ukfMatrix = SquareRootUnscentedKalmanFilter<Matrix<T, 3, 1>>(1,2,1);
    ASSERT_TRUE(ukfMatrix.S.isIdentity()); // S should be identity
}

TEST(SquareRootUnscentedKalmanFilter, computeSigmaPoints) {
    T alpha = 1, beta = 2, kappa = 1;
    
    auto ukf = SquareRootUnscentedKalmanFilter<Vector<T, 3>>(alpha,beta,kappa);
    
    // Init variables
    ukf.gamma = 2;
    ukf.x << 1.f, 2.f, 3.f;
    ukf.S.setIdentity();
    
    ukf.computeSigmaPoints();
    
    ASSERT_FLOAT_EQ(1, ukf.sigmaStatePoints(0,0));
    ASSERT_FLOAT_EQ(2, ukf.sigmaStatePoints(1,0));
    ASSERT_FLOAT_EQ(3, ukf.sigmaStatePoints(2,0));
    
    // Center block diagonal
    ASSERT_FLOAT_EQ(3, ukf.sigmaStatePoints(0,1));
    ASSERT_FLOAT_EQ(4, ukf.sigmaStatePoints(1,2));
    ASSERT_FLOAT_EQ(5, ukf.sigmaStatePoints(2,3));
    
    // Right block diagonal
    ASSERT_FLOAT_EQ(-1, ukf.sigmaStatePoints(0,4));
    ASSERT_FLOAT_EQ(0,  ukf.sigmaStatePoints(1,5));
    ASSERT_FLOAT_EQ(1,  ukf.sigmaStatePoints(2,6));
}

TEST(SquareRootUnscentedKalmanFilter, computeCovarianceSquareRootFromSigmaPoints) {
    T alpha = 1, beta = 2, kappa = 1;
    
    auto ukf = SquareRootUnscentedKalmanFilter<Vector<T, 3>>(alpha,beta,kappa);
    
    // Compute standard (non-square root) covariance from regular UKF
    // and test equality with square-root solution
    
    Cholesky<Matrix<T, 4, 4>> Rsqrt;
    Matrix<T,4,4> _Rsqrt;
    _Rsqrt <<   1, 0, 0, 0,
                0, 2, 0, 0,
                0, 0, 3, 0,
                0, 0, 0, 4;
    
    Rsqrt.setL(_Rsqrt);
    
    Vector<T,4> mean;
    mean << 2, 3, 4, 5;
    
    ukf.sigmaWeights_c << 3, 5, 5, 5, 5, 5, 5;
    
    Matrix<T, 4, 7> sigmaPoints;
    sigmaPoints <<
        1,  2,  1,  1,  0,  1,  1,
        2,  3,  3,  2,  1,  1,  2,
        3,  5,  4,  4,  1,  2,  2,
        3,  4,  5,  2,  0,  4,  3;
    
    Cholesky<Matrix<T,4,4>> S;
    
    ASSERT_TRUE(ukf.computeCovarianceSquareRootFromSigmaPoints(mean, sigmaPoints, Rsqrt, S));
    
    // Compute reference P
    Matrix<T,4,4> P = Matrix<T,4,4>::Zero();
    for( int i = 0; i <= 2*3; ++i )
    {
        Vector<T,4> vec = sigmaPoints.col(i) - mean;
        P += ukf.sigmaWeights_c[i] * vec * vec.transpose();
    }
    P += Rsqrt.reconstructedMatrix();
    
    // Compute squared S
    Matrix<T,4,4> Ssquared = S.reconstructedMatrix().eval();
    
    ASSERT_MATRIX_NEAR(P, Ssquared, 3.5e-5);
}

TEST(SquareRootUnscentedKalmanFilter, computeKalmanGain) {
    T alpha = 1, beta = 2, kappa = 1;
    
    auto ukf = SquareRootUnscentedKalmanFilter<Vector<T, 3>>(alpha,beta,kappa);
    
    ukf.sigmaWeights_c << 3, 5, 5, 5, 5, 5, 5;
    
    Vector<T,2> mean;
    mean << 2, 3;
    
    typename SquareRootUnscentedKalmanFilter<Vector<T, 3>>::template SigmaPoints<decltype(mean)> sigmaPoints;
    sigmaPoints <<
        1,  2,  1,  1,  0,  1,  1,
        2,  3,  3,  2,  1,  1,  2;
    
    Matrix<T,2,2> Ssquared;
    Ssquared <<
        35, 34,
        34, 46;
        
    Cholesky<Matrix<T,2,2>> S;
    S.compute(Ssquared);
    ASSERT_TRUE(S.info() == Eigen::Success);
    
    // x and sigmaStatePoints
    ukf.x << 3, 5, 7;
    ukf.sigmaStatePoints <<
        3,  5, 3, 3, 1, 3, 3,
        5,  7, 7, 5, 3, 3, 5,
        7, 11, 9, 9, 3, 5, 5;
    
    // Reference value for P
    Matrix<T,3,2> P_xy_ref;
    P_xy_ref <<
        20, 20,
        20, 40,
        40, 60;
    
    // Reference value for K
    Matrix<T, 3, 2> K_ref;
    K_ref <<
        0.528634361233480, 0.044052863436123,
       -0.969162995594713, 1.585903083700440,
       -0.440528634361233, 1.629955947136563;
    
    // Kalman gain to be computed
    Matrix<T, 3, 2> K;
    ukf.computeKalmanGain(mean, sigmaPoints, S, K);
    
    ASSERT_MATRIX_NEAR(K_ref, K, 1e-6);
}

TEST(SquareRootUnscentedKalmanFilter, updateStateCovariance) {
    T alpha = 1, beta = 2, kappa = 1;
    
    auto ukf = SquareRootUnscentedKalmanFilter<Vector<T, 3>>(alpha,beta,kappa);
    
    // Setup S_y
    Matrix<T,2,2> Ssquared;
    Ssquared <<
        1, 0.1,
        0.1, 0.5;
    Cholesky<Matrix<T,2,2>> S_y;
    S_y.compute(0.1*Ssquared);
    ASSERT_TRUE(S_y.info() == Eigen::Success);
    
    // Setup Kalman Gain
    Matrix<T, 3, 2> K;
    K <<
        0.178408514951850, -0.304105423213381,
       -1.110998479472884,  0.802838317283324,
       -1.246156445345498,  0.063524243960128;
    
    // Setup S
    ukf.S.setIdentity();
    
    // Setup reference value for S (computed from regular UKF formula)
    Matrix<T,3,3> P_ref = ukf.S.reconstructedMatrix() - K * S_y.reconstructedMatrix() * K.transpose();

    // Perform update
    bool success = ukf.updateStateCovariance<Vector<T,2>>(K, S_y);
    ASSERT_TRUE(success);

    Matrix<T,3,3> P = ukf.S.reconstructedMatrix();
    
    ASSERT_MATRIX_NEAR(P_ref, P, 1e-6);
}
