#ifndef KALMAN_KALMANFILTERBASE_HPP_
#define KALMAN_KALMANFILTERBASE_HPP_

#include "Matrix.hpp"
#include "Types.hpp"

namespace Kalman {
    
    /**
     * @brief Abstract base class for all Kalman Filters
     * 
     * @param StateType The vector-type of the system state (usually some type derived from Kalman::Vector)
     */
    template<class StateType>
    class KalmanFilterBase
    {
    public:
        static_assert( StateType::length >  0, "State vector must contain at least 1 element");
        
        //! Numeric scalar type
        typedef typename StateType::Scalar T;
        
        //! Type of the state vector
        typedef StateType State;
        
    protected:
        //! Estimated state
        State x;
        
    public:
        /**
         * Get current state estimate
         */
        const State& getState() const
        {
            return x;
        }
        
        /**
         * @brief Initialize state
         * @param initialState The initial state of the system
         */
        void init(const State& initialState)
        {
            x = initialState;
        }
    protected:
        /**
         * @brief Protected constructor
         */
        KalmanFilterBase()
        {
        }
    };
}

#endif
