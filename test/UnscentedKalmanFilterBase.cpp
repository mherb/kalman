#include "TestHelper.h"
#include "models/Quadratic.hpp"

#define private public
#define protected public

#include <kalman/UnscentedKalmanFilterBase.hpp>

using namespace Kalman;

template<typename T>
using Vector3 = Vector<T, 3>;

template<typename T>
using Vector2 = Vector<T, 2>;

template<class StateType>
class ConcreteUKF : public UnscentedKalmanFilterBase<StateType>
{
public:
    typedef UnscentedKalmanFilterBase<StateType> Base;
    using typename Base::T;
    
    ConcreteUKF(T alpha, T beta, T kappa) : Base(alpha, beta, kappa) {}
};

typedef float T;

TEST(UnscentedKalmanFilterBase, init) {
    auto ukf = ConcreteUKF<Vector<T, 3>>(1,2,1);
    
    // x should be zero
    ASSERT_FLOAT_EQ(0, ukf.x[0]);
    ASSERT_FLOAT_EQ(0, ukf.x[1]);
    ASSERT_FLOAT_EQ(0, ukf.x[2]);
}

TEST(UnscentedKalmanFilterBase, computeSigmaWeights) {
    T alpha = 1, beta = 2, kappa = 1;
    
    auto ukf = ConcreteUKF<Vector<T, 3>>(alpha,beta,kappa);
    ukf.computeWeights();
    
    ASSERT_FLOAT_EQ(2, ukf.gamma);
    ASSERT_FLOAT_EQ(1, ukf.lambda);
    
    ASSERT_FLOAT_EQ(0.25, ukf.sigmaWeights_m[0]);
    ASSERT_FLOAT_EQ(2.25, ukf.sigmaWeights_c[0]);
    
    for(size_t i = 1; i < 7; i++)
    {
        ASSERT_FLOAT_EQ(0.125, ukf.sigmaWeights_c[i]);
        ASSERT_FLOAT_EQ(0.125, ukf.sigmaWeights_m[i]);
    }
}

TEST(UnscentedKalmanFilterBase, computeSigmaPointTransition) {
    T alpha = 1, beta = 2, kappa = 1;
    
    auto ukf = ConcreteUKF<Vector<T, 3>>(alpha,beta,kappa);

    decltype(ukf.sigmaStatePoints) initial;
    initial <<
            1,  2,  3,  4,  5,  6,  7,
            8,  9, 10, 11, 12, 13, 14,
            15, 16, 17, 18, 19, 20, 21;
    
    // Control vector
    Vector<T,3> u;
    u << 1, 2, 3;
    
    // Compute reference result
    Matrix<T,3,7> ref = (initial.cwiseProduct(initial).colwise() + u).eval();

    // LEGACY System Model
    // Init variables
    ukf.sigmaStatePoints = initial;
    // Create Model
    auto legacyModel = Kalman::Test::Models::QuadraticLinearizedSystemModel<Vector<T, 3>>();
    // Compute transition
    ukf.computeSigmaPointTransition(legacyModel, u);
    ASSERT_MATRIX_NEAR(ref, ukf.sigmaStatePoints, 1e-10);
    // END LEGACY

    // TEMPLATE System Model
    // Init variables
    ukf.sigmaStatePoints = initial;
    // Create Model
    auto templateModel = Kalman::Test::Models::QuadraticTemplateSystemModel<Vector3>();
    // Compute transition
    ukf.computeSigmaPointTransition(templateModel, u);
    ASSERT_MATRIX_NEAR(ref, ukf.sigmaStatePoints, 1e-10);
    // END TEMPLATE
}

TEST(UnscentedKalmanFilterBase, computeSigmaPointMeasurements) {
    T alpha = 1, beta = 2, kappa = 1;
    
    auto ukf = ConcreteUKF<Vector<T, 3>>(alpha,beta,kappa);

    // Init variables
    ukf.sigmaStatePoints <<
         1,  2,  3,  4,  5,  6,  7,
         8,  9, 10, 11, 12, 13, 14,
        15, 16, 17, 18, 19, 20, 21;
    
    // Compute Reference result
    Matrix<T,2,7> tmp = ukf.sigmaStatePoints.template topRows<2>();
    Matrix<T,2,7> ref = tmp.cwiseProduct(tmp).eval();
    
    typename ConcreteUKF<Vector<T,3>>::template SigmaPoints<Vector<T,2>> points;

    // LEGACY System Model
    // Init model
    auto legacyModel = Kalman::Test::Models::QuadraticLinearizedMeasurementModel<Vector<T, 3>, Vector<T, 2>>();
    // Compute points
    ukf.computeSigmaPointMeasurements(legacyModel, points);
    // Compare ref and result
    ASSERT_MATRIX_NEAR(ref, points, 1e-10);
    // END LEGACY

    // TEMPLATE System Model
    // Init model
    auto templateModel = Kalman::Test::Models::QuadraticTemplateMeasurementModel<Vector3, Vector2>();
    // Compute points
    points = decltype(points)::Zero();
    ukf.computeSigmaPointMeasurements<Vector2<T>, decltype(templateModel)>(templateModel, points);
    // Compare ref and result
    ASSERT_MATRIX_NEAR(ref, points, 1e-10);
    // END TEMPLATE
}

TEST(UnscentedKalmanFilterBase, computePredictionFromSigmaPoints) {
    T alpha = 1, beta = 2, kappa = 1;
    
    auto ukf = ConcreteUKF<Vector<T, 3>>(alpha,beta,kappa);
    
    // Init variables
    ukf.sigmaWeights_m << 7, 6, 5, 4, 3, 2, 1;
    
    typename ConcreteUKF<Vector<T,3>>::SigmaPoints<Vector<T,4>> points;
    points <<
         1,  2,  3,  4,  5,  6,  7,
         10,  20,  30,  40,  50,  60,  70,
         100,  200,  300,  400,  500,  600,  700,
         1000,  2000,  3000,  4000,  5000,  6000,  7000;
    
    Vector<T,4> x = ukf.computePredictionFromSigmaPoints<Vector<T,4>>( points );
    
    ASSERT_FLOAT_EQ(84,    x[0]);
    ASSERT_FLOAT_EQ(840,   x[1]);
    ASSERT_FLOAT_EQ(8400,  x[2]);
    ASSERT_FLOAT_EQ(84000, x[3]);
}
