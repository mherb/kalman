#include "TestHelper.h"

#define private public
#define protected public

#include <kalman/Model.hpp>
#include "models/Quadratic.hpp"

using namespace Kalman;

template<typename T>
using Vector3 = Vector<T, 3>;

typedef Vector3<float> Vector3f;
typedef Vector3<double> Vector3d;

TEST(Model, getCovariance) {
    // These tests are primarly to check for proper inference and compilation of the templates
    {
        // Get covariance with differing scalar types (implicit casting)
        typedef Kalman::Test::Models::QuadraticTemplateJacobianSystemModel<float, Vector3> System;
        System system;
        Covariance<double, Vector3d> cov = Model::getCovariance<System, double, Vector3d>(system);
        ASSERT_MATRIX_DOUBLE_EQ(decltype(cov)::Identity(), cov);
    }
    {
        // Get covariance from templated function
        typedef Kalman::Test::Models::QuadraticTemplateSystemModel<Vector3> System;
        System system;
        Covariance<double, Vector3d> cov = Model::getCovariance<System, double, Vector3d>(system);
        ASSERT_MATRIX_DOUBLE_EQ(decltype(cov)::Identity(), cov);
    }
}

TEST(Model, hasJacobian) {
    {
        typedef Kalman::Test::Models::QuadraticTemplateJacobianSystemModel<float, Vector3> System;
        bool hasJacobian = Model::hasJacobian<System, Vector3f /* state */, Vector3f /* control */>::value;
        ASSERT_TRUE(hasJacobian);
    }
    {
        typedef Kalman::Test::Models::QuadraticTemplateSystemModel<Vector3> System;
        bool hasJacobian = Model::hasJacobian<System, Vector3f /* state */, Vector3f /* control */>::value;
        ASSERT_FALSE(hasJacobian);
    }
}

TEST(Model, isTemplateModel) {
    {
        typedef Kalman::Test::Models::QuadraticTemplateJacobianSystemModel<float, Vector3> System;
        bool isTemplate = Model::isTemplateModel<System, Vector3f /* state */, Vector3f /* control */>::value;
        ASSERT_TRUE(isTemplate);
    }
    {
        typedef Kalman::Test::Models::QuadraticTemplateJacobianMeasurementModel<float, Vector3, Vector3> System;
        bool isTemplate = Model::isTemplateModel<System, Vector3f /* state */, Vector3f /* measurement */>::value;
        ASSERT_TRUE(isTemplate);
    }
    {
        typedef Kalman::Test::Models::QuadraticTemplateSystemModel<Vector3> System;
        bool isTemplate = Model::isTemplateModel<System, Vector3f /* state */, Vector3f /* control */>::value;
        ASSERT_TRUE(isTemplate);
    }
    {
        typedef Kalman::Test::Models::QuadraticTemplateMeasurementModel<Vector3, Vector3> System;
        bool isTemplate = Model::isTemplateModel<System, Vector3f /* state */, Vector3f /* measurement */>::value;
        ASSERT_TRUE(isTemplate);
    }
    {
        typedef Kalman::Test::Models::QuadraticLinearizedSystemModel<Vector3f> System;
        bool isTemplate = Model::isTemplateModel<System, Vector3f /* state */, Vector3f /* control */>::value;
        ASSERT_FALSE(isTemplate);
    }
    {
        typedef Kalman::Test::Models::QuadraticLinearizedMeasurementModel<Vector3f> System;
        bool isTemplate = Model::isTemplateModel<System, Vector3f /* state */, Vector3f /* measurement */>::value;
        ASSERT_FALSE(isTemplate);
    }
    {
        typedef Kalman::Test::Models::QuadraticCombinedJacobianSystemModel<float, Vector3> System;
        bool isTemplate = Model::isTemplateModel<System, Vector3f /* state */, Vector3f /* control */>::value;
        ASSERT_TRUE(isTemplate);
    }
    {
        typedef Kalman::Test::Models::QuadraticCombinedJacobianMeasurementModel<float, Vector3, Vector3> System;
        bool isTemplate = Model::isTemplateModel<System, Vector3f /* state */, Vector3f /* measurement */>::value;
        ASSERT_TRUE(isTemplate);
    }
}

TEST(Model, isLegacyModel) {
    {
        typedef Kalman::Test::Models::QuadraticTemplateJacobianSystemModel<float, Vector3> System;
        bool legacy = Model::isLegacyModel<System, Vector3f /* state */, Vector3f /* control */>::value;
        ASSERT_FALSE(legacy);
    }
    {
        typedef Kalman::Test::Models::QuadraticTemplateJacobianMeasurementModel<float, Vector3, Vector3> System;
        bool legacy = Model::isLegacyModel<System, Vector3f /* state */>::value;
        ASSERT_FALSE(legacy);
    }
    {
        typedef Kalman::Test::Models::QuadraticTemplateSystemModel<Vector3> System;
        bool legacy = Model::isLegacyModel<System, Vector3f /* state */, Vector3f /* control */>::value;
        ASSERT_FALSE(legacy);
    }
    {
        typedef Kalman::Test::Models::QuadraticTemplateMeasurementModel<Vector3, Vector3> System;
        bool legacy = Model::isLegacyModel<System, Vector3f /* state */>::value;
        ASSERT_FALSE(legacy);
    }
    {
        typedef Kalman::Test::Models::QuadraticLinearizedSystemModel<Vector3f> System;
        bool legacy = Model::isLegacyModel<System, Vector3f /* state */, Vector3f /* control */>::value;
        ASSERT_TRUE(legacy);
    }
    {
        typedef Kalman::Test::Models::QuadraticLinearizedMeasurementModel<Vector3f> System;
        bool legacy = Model::isLegacyModel<System, Vector3f /* state */>::value;
        ASSERT_TRUE(legacy);
    }
    {
        typedef Kalman::Test::Models::QuadraticCombinedJacobianSystemModel<float, Vector3> System;
        bool legacy = Model::isLegacyModel<System, Vector3f /* state */, Vector3f /* control */>::value;
        ASSERT_TRUE(legacy);
    }
    {
        typedef Kalman::Test::Models::QuadraticCombinedJacobianMeasurementModel<float, Vector3, Vector3> System;
        bool legacy = Model::isLegacyModel<System, Vector3f /* state */>::value;
        ASSERT_TRUE(legacy);
    }
}
