#include "TestHelper.h"

#define private public
#define protected public

#include <kalman/UnscentedKalmanFilter.hpp>

using namespace Kalman;

typedef float T;

TEST(UnscentedKalmanFilter, init) {
    auto ukf = UnscentedKalmanFilter<Vector<T, 3>>(1,2,1);
    ASSERT_TRUE(ukf.P.isIdentity()); // P should be identity

    // Same as above, but with general matrix type instead of vector
    auto ukfMatrix = UnscentedKalmanFilter<Matrix<T, 3, 1>>(1,2,1);
    ASSERT_TRUE(ukfMatrix.P.isIdentity()); // P should be identity
}

TEST(UnscentedKalmanFilter, computeSigmaPoints) {
    T alpha = 1, beta = 2, kappa = 1;
    
    auto ukf = UnscentedKalmanFilter<Vector<T, 3>>(alpha,beta,kappa);
    
    // Init variables
    ukf.gamma = 2;
    ukf.x << 1.f, 2.f, 3.f;
    ukf.P.setIdentity();
    
    ASSERT_TRUE(ukf.computeSigmaPoints());
    
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

TEST(UnscentedKalmanFilter, computeCovarianceSquareRootFromSigmaPoints) {
    T alpha = 1, beta = 2, kappa = 1;
    
    auto ukf = UnscentedKalmanFilter<Vector<T, 3>>(alpha,beta,kappa);
    
    Matrix<T, 4, 4> R;
    R <<    1, 0, 0, 0,
            0, 4, 0, 0,
            0, 0, 9, 0,
            0, 0, 0, 16;
    
    Vector<T,4> mean;
    mean << 2, 3, 4, 5;
    
    ukf.sigmaWeights_c << 3, 5, 5, 5, 5, 5, 5;
    
    Matrix<T, 4, 7> sigmaPoints;
    sigmaPoints <<
        1,  2,  1,  1,  0,  1,  1,
        2,  3,  3,  2,  1,  1,  2,
        3,  5,  4,  4,  1,  2,  2,
        3,  4,  5,  2,  0,  4,  3;
    
    Matrix<T,4,4> P;
    
    ASSERT_TRUE(ukf.computeCovarianceFromSigmaPoints(mean, sigmaPoints, R, P));
    
    Matrix<T,4,4> P_ref;
    P_ref <<
        44, 43,  53,  86,
        43, 57,  63,  91,
        53, 63, 102, 106,
        86, 91, 106, 228;
    
    ASSERT_MATRIX_NEAR(P_ref, P, 1e-6);
}

TEST(UnscentedKalmanFilter, computeKalmanGain) {
    T alpha = 1, beta = 2, kappa = 1;
    
    auto ukf = UnscentedKalmanFilter<Vector<T, 3>>(alpha,beta,kappa);
    
    ukf.sigmaWeights_c << 3, 5, 5, 5, 5, 5, 5;
    
    Vector<T,2> mean;
    mean << 2, 3;
    
    typename UnscentedKalmanFilter<Vector<T, 3>>::template SigmaPoints<decltype(mean)> sigmaPoints;
    sigmaPoints <<
        1,  2,  1,  1,  0,  1,  1,
        2,  3,  3,  2,  1,  1,  2;
    
    Matrix<T,2,2> P_yy;
    P_yy <<
        35, 34,
        34, 46;
    
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
    ukf.computeKalmanGain(mean, sigmaPoints, P_yy, K);
    
    ASSERT_MATRIX_NEAR(K_ref, K, 1e-6);
}

TEST(UnscentedKalmanFilter, updateStateCovariance) {
    T alpha = 1, beta = 2, kappa = 1;
    
    auto ukf = UnscentedKalmanFilter<Vector<T, 3>>(alpha,beta,kappa);
    
    Matrix<T,2,2> P_yy;
    P_yy <<
        1, 0.1,
        0.1, 0.5;
    P_yy *= 0.1;
    
    // Setup Kalman Gain
    Matrix<T, 3, 2> K;
    K <<
        0.178408514951850, -0.304105423213381,
       -1.110998479472884,  0.802838317283324,
       -1.246156445345498,  0.063524243960128;
    
    // Setup P
    ukf.P.setIdentity();
    
    // Setup P_ref
    Matrix<T, 3, 3> P_ref;
    P_ref <<
        0.993278134696764,   0.027217594648895,   0.019295433443564,
        0.027217594648895,   0.862179812671265,  -0.130287401631777,
        0.019295433443564,  -0.130287401631777,   0.846090867814784;
    
    // Perform update
    bool success = ukf.updateStateCovariance<Vector<T,2>>(K, P_yy);
    ASSERT_TRUE(success);
    
    ASSERT_MATRIX_NEAR(P_ref, ukf.P, 1e-6);
}
