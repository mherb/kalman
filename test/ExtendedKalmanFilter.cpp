#include "TestHelper.h"

#define private public
#define protected public

#include <kalman/ExtendedKalmanFilter.hpp>

using namespace Kalman;

typedef float T;

TEST(ExtendedKalmanFilter, init) {
    ExtendedKalmanFilter<Vector<T, 3>> ekf;
    ASSERT_TRUE(ekf.P.isIdentity()); // P should be identity

    // Same as above, but with general matrix type instead of vector
    ExtendedKalmanFilter<Matrix<T, 3, 1>> ekfMatrix;
    ASSERT_TRUE(ekfMatrix.P.isIdentity()); // P should be identity
}
