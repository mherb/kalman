#ifndef KALMAN_SQUAREROOTFILTERBASE_HPP_
#define KALMAN_SQUAREROOTFILTERBASE_HPP_

#include "SquareRootBase.hpp"

namespace Kalman {
    
    /**
     * @brief Abstract base class for square root filters
     * 
     * @param StateType The vector-type of the system state (usually some type derived from Kalman::Vector)
     */
    template<class StateType>
    class SquareRootFilterBase : public SquareRootBase<StateType>
    {
    protected:
        //! SquareRoot Base Type
        typedef SquareRootBase<StateType> Base;
        
        //! Covariance Square Root
        using Base::S;
    };
}

#endif
