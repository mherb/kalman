#ifndef KALMAN_LINEARIZEDMEASUREMENTMODEL_HPP_
#define KALMAN_LINEARIZEDMEASUREMENTMODEL_HPP_

#include "MeasurementModel.hpp"

namespace Kalman {
    template<class StateType>
    class ExtendedKalmanFilter;
    template<class StateType>
    class SquareRootExtendedKalmanFilter;
    
    /**
     * @brief Abstract base class of all linearized (first order taylor expansion) measurement models
     *
     * @param StateType The vector-type of the system state (usually some type derived from Kalman::Vector)
     * @param MeasurementType The vector-type of the measurement (usually some type derived from Kalman::Vector)
     * @param CovarianceBase The class template used for covariance storage (must be either StandardBase or SquareRootBase)
     */
    template<class StateType, class MeasurementType, template<class> class CovarianceBase = StandardBase>
    class LinearizedMeasurementModel : public MeasurementModel<StateType, MeasurementType, CovarianceBase>
    {
        friend class ExtendedKalmanFilter<StateType>;
        friend class SquareRootExtendedKalmanFilter<StateType>;
    public:
        //! Measurement model base
        typedef MeasurementModel<StateType, MeasurementType, CovarianceBase> Base;
        
        //! System state type
        using typename Base::State;
        
        //! Measurement vector type
        using typename Base::Measurement;
        
    protected:
        //! Measurement model jacobian
        Jacobian<Measurement, State> H;
        //! Measurement model noise jacobian
        Jacobian<Measurement, Measurement> V;
        
        /**
         * Callback function for state-dependent update of Jacobi-matrices H and V before each update step
         */
        virtual void updateJacobians( const State& x )
        {
            // No update by default
        }
    protected:
        LinearizedMeasurementModel()
        {
            H.setIdentity();
            V.setIdentity();
        }
        ~LinearizedMeasurementModel() {}
    };
}

#endif
