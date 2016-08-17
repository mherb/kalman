#include "TestHelper.h"

#define private public
#define protected public

#include <kalman/SquareRootBase.hpp>

using namespace Kalman;

typedef float T;
typedef Kalman::Vector<T, 3> Vec3f;

TEST(SquareRootBase, setCovariance) {
    Kalman::Covariance<Vec3f> cov;
    cov <<
    1, 2, 4,
    2, 13, 23,
    4, 23, 77;

    SquareRootBase<Vec3f> S;
    S.setCovariance(cov);

    ASSERT_MATRIX_NEAR(cov, S.S.reconstructedMatrix(), 1e-5);
}
