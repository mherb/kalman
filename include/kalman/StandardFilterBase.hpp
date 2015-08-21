#ifndef KALMAN_STANDARDFILTERBASE_HPP_
#define KALMAN_STANDARDFILTERBASE_HPP_

#include "StandardBase.hpp"

namespace Kalman {
    
    /**
     * @brief Abstract base class for standard (non-square root) filters
     * 
     * @param StateType The vector-type of the system state (usually some type derived from Kalman::Vector)
     */
    template<class StateType>
    class StandardFilterBase : public StandardBase<StateType>
    {
    protected:
        //! Standard Base Type
        typedef StandardBase<StateType> Base;
        
        //! Covariance matrix
        using Base::P;
    };
}

#endif
