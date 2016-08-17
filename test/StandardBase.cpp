#include "TestHelper.h"

#define private public
#define protected public

#include <kalman/StandardBase.hpp>

using namespace Kalman;

typedef float T;
typedef Kalman::Vector<T, 3> Vec3f;

TEST(StandardBase, setCovarianceSquareRoot) {
    Kalman::Covariance<Vec3f> sqrRoot;
    sqrRoot <<
    1, 0, 0,
    2, 3, 0,
    4, 5, 6;

    Kalman::Covariance<Vec3f> cov;
    cov <<
    1, 2, 4,
    2, 13, 23,
    4, 23, 77;

    StandardBase<Vec3f> S;
    S.setCovarianceSquareRoot(sqrRoot);

    ASSERT_MATRIX_NEAR(cov, S.P, 1e-5);
}
