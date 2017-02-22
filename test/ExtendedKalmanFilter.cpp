#include "TestHelper.h"

#define private public
#define protected public

#include <kalman/ExtendedKalmanFilter.hpp>
#include "models/Quadratic.hpp"

using namespace Kalman;

typedef float T;
template<typename _T>
using Vector3 = Vector<_T, 3>;
template<typename _T>
using Vector2 = Vector<_T, 2>;

TEST(ExtendedKalmanFilter, init) {
    ExtendedKalmanFilter<Vector3<T>> ekf;
    ASSERT_TRUE(ekf.P.isIdentity()); // P should be identity
}

TEST(ExtendedKalmanFilter, predict) {
    ExtendedKalmanFilter<Vector3<T>> ekf;

    Vector3<T> initial;
    initial << 1, 2, 3;

    Vector3<T> u;
    u << 0, 1, 2;

    // LEGACY Linearized System Model
    Kalman::Test::Models::QuadraticLinearizedSystemModel<Vector3<T>> sysLinearized;

    ekf.init(initial);
    ekf.P = decltype(ekf.P)::Identity();
    ekf.predict(sysLinearized, u);

    ASSERT_FLOAT_EQ(1,  ekf.x.x());
    ASSERT_FLOAT_EQ(5,  ekf.x.y());
    ASSERT_FLOAT_EQ(11, ekf.x.z());
    // END LEGACY

    // AUTO-DIFF
#ifdef KALMAN_EKF_SUPPORTS_AUTODIFF
    Kalman::Test::Models::QuadraticTemplateSystemModel<Vector3> sysAutoDiff;

    ekf.init(initial);
    ekf.P = decltype(ekf.P)::Identity();
    ekf.predict(sysAutoDiff, u);

    ASSERT_FLOAT_EQ(1,  ekf.x.x());
    ASSERT_FLOAT_EQ(5,  ekf.x.y());
    ASSERT_FLOAT_EQ(11, ekf.x.z());
    // END AUTO-DIFF
#else
    std::cerr << "ExtendedKalmanFilter: AutoDiff support is disabled because it is not supported by the compiler!" << std::endl;
#endif

    // EXPLICIT JACOBIAN
    Kalman::Test::Models::QuadraticTemplateJacobianSystemModel<T, Vector3> sysJacobian;

    ekf.init(initial);
    ekf.P = decltype(ekf.P)::Identity();
    ekf.predict(sysJacobian, u);

    ASSERT_FLOAT_EQ(1,  ekf.x.x());
    ASSERT_FLOAT_EQ(5,  ekf.x.y());
    ASSERT_FLOAT_EQ(11, ekf.x.z());
    // END EXPLICIT JACOBIAN
}

TEST(ExtendedKalmanFilter, update) {
    ExtendedKalmanFilter<Vector3<T>> ekf;

    Vector3<T> initial;
    initial << 1, 2, 3;
    ekf.init(initial);

    Vector2<T> z;
    z << 0.5, 4;

    Vector2<T> innovation = z - initial.cwiseProduct(initial).template head<2>();
    Vector3<T> expected = initial;
    expected.x() += innovation.x() * 2*initial.x()/(1+4*initial.x()*initial.x());
    expected.y() += innovation.y() * 2*initial.y()/(1+4*initial.y()*initial.y());

    // LEGACY Linearized Measurement Model
    Kalman::Test::Models::QuadraticLinearizedMeasurementModel<Vector3<T>, Vector2<T>> measLinearized;

    ekf.init(initial);
    ekf.P = decltype(ekf.P)::Identity();
    ekf.update(measLinearized, z);

    ASSERT_FLOAT_EQ(expected.x(),  ekf.x.x());
    ASSERT_FLOAT_EQ(expected.y(),  ekf.x.y());
    ASSERT_FLOAT_EQ(initial.z(), ekf.x.z());
    // END LEGACY

    // AUTO-DIFF
#ifdef KALMAN_EKF_SUPPORTS_AUTODIFF
    Kalman::Test::Models::QuadraticTemplateMeasurementModel<Vector3, Vector2> measAutoDiff;

    ekf.init(initial);
    ekf.P = decltype(ekf.P)::Identity();
    ekf.update(measAutoDiff, z);

    ASSERT_FLOAT_EQ(expected.x(),  ekf.x.x());
    ASSERT_FLOAT_EQ(expected.y(),  ekf.x.y());
    ASSERT_FLOAT_EQ(initial.z(), ekf.x.z());
    // END AUTO-DIFF
#else
    std::cerr << "ExtendedKalmanFilter: AutoDiff support is disabled because it is not supported by the compiler!" << std::endl;
#endif

    // EXPLICIT JACOBIAN
    Kalman::Test::Models::QuadraticTemplateJacobianMeasurementModel<T, Vector3, Vector2> measJacobian;

    ekf.init(initial);
    ekf.P = decltype(ekf.P)::Identity();
    ekf.update(measJacobian, z);

    ASSERT_FLOAT_EQ(expected.x(),  ekf.x.x());
    ASSERT_FLOAT_EQ(expected.y(),  ekf.x.y());
    ASSERT_FLOAT_EQ(initial.z(), ekf.x.z());
    // END EXPLICIT JACOBIAN
}
