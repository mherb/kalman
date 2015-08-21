#ifndef KALMAN_STANDARDBASE_HPP_
#define KALMAN_STANDARDBASE_HPP_

#include "Types.hpp"

namespace Kalman {
    
    /**
     * @brief Abstract base class for standard (non-square root) filters and models
     * 
     * @param StateType The vector-type of the system state (usually some type derived from Kalman::Vector)
     */
    template<class StateType>
    class StandardBase
    {
    protected:
        //! Covariance
        Covariance<StateType> P;
        
    public:
        /**
         * Get covariance
         */
        const Covariance<StateType>& getCovariance() const
        {
            return P;
        }
        
        /**
         * Get covariance (as square root)
         */
        CovarianceSquareRoot<StateType> getCovarianceSquareRoot() const
        {
            return CovarianceSquareRoot<StateType>(P);
        }
        
        /**
         * Set Covariance
         */
        bool setCovariance(const Covariance<StateType>& covariance)
        {
            P = covariance;
            return true;
        }
        
    protected:
        StandardBase()
        {
            P.setIdentity();
        }
    };
}

#endif
