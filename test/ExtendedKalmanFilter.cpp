#include "TestHelper.h"

#define private public
#define protected public

#include <kalman/ExtendedKalmanFilter.hpp>

using namespace Kalman;

typedef float T;

TEST(ExtendedKalmanFilter, init) {
    ExtendedKalmanFilter<Vector<T, 3>> ekf;
    
    // P should be identity
    ASSERT_TRUE(ekf.P.isIdentity());
}
