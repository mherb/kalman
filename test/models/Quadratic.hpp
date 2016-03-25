#ifndef KALMAN_TEST_MODELS_QUADRATIC_HPP_
#define KALMAN_TEST_MODELS_QUADRATIC_HPP_

#include <kalman/LinearizedSystemModel.hpp>
#include <kalman/LinearizedMeasurementModel.hpp>

namespace Kalman
{
namespace Test
{
namespace Models
{

template<class StateType, template<typename> class CovarianceBase = Kalman::StandardBase>
class QuadraticLinearizedSystemModel : public Kalman::LinearizedSystemModel<StateType, StateType, CovarianceBase>
{
public:
    typedef SystemModel<StateType, StateType> Base;
    using typename Base::State;
    using typename Base::Control;

    State f(const State& x, const Control& u) const
    {
        // return x.^2 + u
        return x.cwiseProduct(x) + u;
    }

    void updateJacobians(const State& x, const Control& u)
    {
        this->F = (2. * x.asDiagonal()).toDenseMatrix();
        (void)u;
    }
};

template<class StateType, class MeasurementType = StateType, template<typename> class CovarianceBase = Kalman::StandardBase>
class QuadraticLinearizedMeasurementModel : public Kalman::LinearizedMeasurementModel<StateType, MeasurementType, CovarianceBase>
{
public:
    typedef MeasurementModel<StateType, MeasurementType> Base;
    using typename Base::State;
    using typename Base::Measurement;

    static_assert(static_cast<decltype(Kalman::Dynamic)>(MeasurementType::RowsAtCompileTime) <= static_cast<decltype(Kalman::Dynamic)>(StateType::RowsAtCompileTime),
                  "Measurement length must be less than or equal to State length");

    Measurement h(const State& x) const
    {
        // return x.^2
        return x.cwiseProduct(x).template head<Measurement::RowsAtCompileTime>();
    }

    void updateJacobians(const State& x)
    {
        this->H = (2. * x.asDiagonal()).toDenseMatrix().template topLeftCorner<Measurement::RowsAtCompileTime, State::RowsAtCompileTime>();
    }
};

} // namespace Models
} // namespace Test
} // namespace Kalman

#endif
