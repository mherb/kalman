#ifndef KALMAN_LINEARIZEDSYSTEMMODEL_HPP_
#define KALMAN_LINEARIZEDSYSTEMMODEL_HPP_

#include "SystemModel.hpp"

namespace Kalman {
    template<class StateType>
    class ExtendedKalmanFilter;
    template<class StateType>
    class SquareRootExtendedKalmanFilter;
    
    /**
     * @brief Abstract base class of all linearized (first order taylor expansion) system models
     *
     * @param StateType The vector-type of the system state (usually some type derived from Kalman::Vector)
     * @param ControlType The vector-type of the control input (usually some type derived from Kalman::Vector)
     * @param CovarianceBase The class template used for covariance storage (must be either StandardBase or SquareRootBase)
     */
    template<class StateType, class ControlType = Vector<typename StateType::Scalar, 0>, template<class> class CovarianceBase = StandardBase >
    class LinearizedSystemModel : public SystemModel<StateType, ControlType, CovarianceBase>
    {
        friend class ExtendedKalmanFilter<StateType>;
        friend class SquareRootExtendedKalmanFilter<StateType>;
    public:
        //! System model base
        typedef SystemModel<StateType, ControlType, CovarianceBase> Base;
        
        //! System state type
        using typename Base::State;
        
        //! System control input type
        using typename Base::Control;
        
    protected:
        //! System model jacobian
        Jacobian<State, State> F;
        //! System model noise jacobian
        Jacobian<State, State> W;
        
        /**
         * Callback function for state-dependent update of Jacobi-matrices F and W before each update step
         */
        virtual void updateJacobians( const State& x, const Control& u )
        {
            // No update by default
        }
    protected:
        LinearizedSystemModel()
        {
            F.setIdentity();
            W.setIdentity();
        }
        ~LinearizedSystemModel() {}
    };
}

#endif
