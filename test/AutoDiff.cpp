#include "TestHelper.h"

#define private public
#define protected public

#include <kalman/AutoDiff.hpp>
#include <kalman/Matrix.hpp>

using namespace Kalman;

template<typename T>
using Vector3 = Vector<T, 3>;


TEST(AutoDiff, evaluate) {
    using ADJ = AutoDiffJacobian<double, Vector3, Vector3>;

    ADJ::Function func = [](const ADJ::ActiveInput& in, ADJ::ActiveOutput& out) {
        out[0] = in[0];
        out[1] = ADJ::ActiveScalar(2)*in[1]*in[2];
        out[2] = sin(3*in[0]+in[2]);
    };

    ADJ diff1(func);

    ADJ::Jacobian expectedJacobian, actualJacobian;
    Vector3<double> input, expectedResult, actualResult;

    input << 1, 2, 3;

    diff1(input, actualResult, actualJacobian);

    expectedResult << input[0], 2*input[1]*input[2], sin(3*input[0]+input[2]);
    expectedJacobian << 1, 0, 0,
                        0, 2*input[2], 2*input[1],
                        cos(3*input[0]+input[2])*3, 0, cos(3*input[0]+input[2]);

    ASSERT_MATRIX_DOUBLE_EQ(expectedJacobian, actualJacobian);
    ASSERT_MATRIX_DOUBLE_EQ(expectedResult, actualResult);
}
