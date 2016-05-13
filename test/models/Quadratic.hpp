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

template<template<typename> class StateType>
class QuadraticTemplateSystemModel {
public:

    template<typename T, typename _T>
    void predict(const StateType<T>& x, const StateType<_T>& u,
                 StateType<T>& prediction) {
        // return x.^2 + u
        prediction = x.cwiseProduct(x) + u.template cast<T>();
    }

    template<typename T>
    Covariance <T, StateType<T> > getCovariance()
    {
        return Covariance <T, StateType<T> >::Identity();
    }
};

template<typename T, template<typename> class StateType>
class QuadraticTemplateJacobianSystemModel {
public:

    void predict(const StateType<T>& x, const StateType<T>& u,
                 StateType<T>& prediction)
    {
        // return x.^2 + u
        prediction = x.cwiseProduct(x) + u;
    }

    Jacobian<T, StateType<T>> getJacobian(const StateType<T>& x,
                                          const StateType<T>& u)
    {
        return 2. * x.asDiagonal();
    }

    Covariance <T, StateType<T> > getCovariance()
    {
        return Covariance<T, StateType<T> >::Identity();
    }
};

template<template<typename> class StateType, template<typename> class MeasurementType>
class QuadraticTemplateMeasurementModel {
public:

    template<typename T>
    void measure(const StateType<T>& x, MeasurementType<T>& measurement) {
        // return x.^2
        measurement = x.cwiseProduct(x).template head<MeasurementType<T>::RowsAtCompileTime>();
    }

    template<typename T>
    Covariance <T, MeasurementType<T> > getCovariance()
    {
        return Covariance<T, MeasurementType<T> >::Identity();
    }
};

template<typename T, template<typename> class StateType, template<typename> class MeasurementType>
class QuadraticTemplateJacobianMeasurementModel {
public:
    typedef StateType<T> State;
    typedef MeasurementType<T> Measurement;

    typedef Kalman::Jacobian<T, State, Measurement> Jacobian;
public:
    void measure(const State& x, Measurement& measurement) {
        // return x.^2
        measurement = x.cwiseProduct(x).template head<MeasurementType<T>::RowsAtCompileTime>();
    }

    Jacobian getJacobian(const StateType<T>& x)
    {
        return (2. * x.asDiagonal()).toDenseMatrix().template topLeftCorner<Measurement::RowsAtCompileTime, State::RowsAtCompileTime>();
    }

    Covariance <T, Measurement> getCovariance()
    {
        return Covariance<T, Measurement>::Identity();
    }
};

template<typename T, template<typename> class StateType>
class QuadraticCombinedJacobianSystemModel :
    public QuadraticLinearizedSystemModel<StateType<T>>,
    public QuadraticTemplateJacobianSystemModel<T, StateType>
{
};

template<typename T, template<typename> class StateType, template<typename> class MeasurementType>
class QuadraticCombinedJacobianMeasurementModel :
    public QuadraticLinearizedMeasurementModel<StateType<T>, MeasurementType<T>>,
    public QuadraticTemplateJacobianMeasurementModel<T, StateType, MeasurementType>
{
};

} // namespace Models
} // namespace Test
} // namespace Kalman

#endif
