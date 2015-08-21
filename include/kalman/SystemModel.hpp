#ifndef KALMAN_SYSTEMMODEL_HPP_
#define KALMAN_SYSTEMMODEL_HPP_

#include <type_traits>

#include "Matrix.hpp"
#include "StandardBase.hpp"

namespace Kalman {
    /**
     * @brief Abstract base class of all system models
     *
     * @param StateType The vector-type of the system state (usually some type derived from Kalman::Vector)
     * @param ControlType The vector-type of the control input (usually some type derived from Kalman::Vector)
     * @param CovarianceBase The class template used for covariance storage (must be either StandardBase or SquareRootBase)
     */
    template<class StateType, class ControlType = Vector<typename StateType::Scalar, 0>, template<class> class CovarianceBase = StandardBase>
    class SystemModel : public CovarianceBase<StateType>
    {
        static_assert( StateType::length   >  0, "State vector must contain at least 1 element");
        static_assert( ControlType::length >= 0, "Control vector must contain at least 0 elements");
        static_assert( std::is_same<typename StateType::Scalar, typename ControlType::Scalar>::value, "State and Control scalar types must be identical" );
    public:
        //! System state type
        typedef StateType State;
        
        //! System control input type
        typedef ControlType Control;
        
    public:
        /**
         * Control Model Function h
         * 
         * Predicts the estimated measurement value given the current state estimate x
         */
        virtual State f(const State& x, const Control& u) const = 0;
        
    protected:
        SystemModel() {}
    };
}

#endif
