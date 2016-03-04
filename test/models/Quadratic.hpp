#ifndef KALMAN_TEST_MODELS_QUADRATIC_HPP_
#define KALMAN_TEST_MODELS_QUADRATIC_HPP_

#include <kalman/SystemModel.hpp>
#include <kalman/MeasurementModel.hpp>

namespace Kalman
{
namespace Test
{
namespace Models
{

template<class StateType>
class QuadraticSystemModel : public SystemModel<StateType, StateType>
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
};

template<class StateType, class MeasurementType = StateType>
class QuadraticMeasurementModel : public MeasurementModel<StateType, MeasurementType>
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
};

} // namespace Models
} // namespace Test
} // namespace Kalman

#endif
