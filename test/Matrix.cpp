#include <gtest/gtest.h>

#define private public
#define protected public

#include <kalman/Matrix.hpp>

using namespace Kalman;

typedef float T;

TEST(Cholesky, solve) {
    Matrix<T, 4, 3> K;
    Matrix<T, 4, 3> P;
    
    // Set S (lower choleksy factor)
    Matrix<T,3,3> S_;
    S_ << 1, 0, 0, 1, 1, 0, 1, 1, 2;
    
    Cholesky<Matrix<T,3,3>> S;
    S.setL(S_);
    
    // Set P
    P << 6, 4, 2, 4, 8, 6, 2, 6, 0, 3, 5, 7;
    
    // Compute K
    K = S.solve(P.transpose()).transpose();

    Matrix<T, 4, 3> P_ = K * S.matrixL().toDenseMatrix() * S.matrixU().toDenseMatrix();
    
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            ASSERT_FLOAT_EQ( P(i,j), P_(i,j) );
        }
    }
}

TEST(QR, compute) {
    Matrix<T, 3, 5> A;
    
    A <<
    8, 7, 1, 1,  4,
    8, 2, 3, 9, 10,
    4, 8, 1, 7,  1;
    
    Eigen::HouseholderQR<Matrix<T,5,3>> qr(A.transpose());
    
    auto R = qr.matrixQR();
    
    ASSERT_FLOAT_EQ(-11.445523142259598f, R(0,0));
    ASSERT_FLOAT_EQ(-11.358152736593492f, R(0,1));
    ASSERT_FLOAT_EQ(-8.737040566610379f, R(0,2));
    ASSERT_FLOAT_EQ(11.357480636664706f, R(1,1));
    ASSERT_FLOAT_EQ(2.180356680396514f, R(1,2));
    ASSERT_FLOAT_EQ(-7.064712795553326f, R(2,2));
}

