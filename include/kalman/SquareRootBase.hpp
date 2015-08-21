#ifndef KALMAN_SQUAREROOTBASE_HPP_
#define KALMAN_SQUAREROOTBASE_HPP_

#include "Types.hpp"

namespace Kalman {
    
    /**
     * @brief Abstract base class for square-root filters and models
     * 
     * @param StateType The vector-type of the system state (usually some type derived from Kalman::Vector)
     */
    template<class StateType>
    class SquareRootBase
    {
    protected:
        //! Covariance Square Root
        CovarianceSquareRoot<StateType> S;
        
    public:
        /**
         * Get covariance (as square root)
         */
        const CovarianceSquareRoot<StateType>& getCovarianceSquareRoot() const
        {
            return S;
        }
        
        /**
         * Get covariance reconstructed from square root
         */
        Covariance<StateType> getCovariance() const
        {
            return S.reconstructedMatrix();
        }
        
        /**
         * Set Covariance
         */
        bool setCovariance(const Covariance<StateType>& covariance)
        {
            S.compute(covariance);
            return (S.info() == Eigen::Success);
        }
        
    protected:
        SquareRootBase()
        {
            S.setIdentity();
        }
    };
}

#endif
