#ifndef KALMAN_EXAMPLES_ROBOT1_POSITIONMEASUREMENTMODEL_HPP_
#define KALMAN_EXAMPLES_ROBOT1_POSITIONMEASUREMENTMODEL_HPP_

#include <kalman/LinearizedMeasurementModel.hpp>

namespace KalmanExamples
{
namespace Robot1
{

/**
 * @brief Measurement vector measuring the robot position
 *
 * @param T Numeric scalar type
 */
template<typename T>
class PositionMeasurement : public Kalman::Vector<T, 2>
{
public:
    KALMAN_VECTOR(PositionMeasurement, T, 2)
    
    //! Distance to landmark 1
    static constexpr size_t D1 = 0;
    
    //! Distance to landmark 2
    static constexpr size_t D2 = 1;
    
    T d1()       const { return (*this)[ D1 ]; }
    T d2()       const { return (*this)[ D2 ]; }
    
    T& d1()      { return (*this)[ D1 ]; }
    T& d2()      { return (*this)[ D2 ]; }
};

/**
 * @brief Measurement model for measuring the position of the robot
 *        using two beacon-landmarks
 *
 * This is the measurement model for measuring the position of the robot.
 * The measurement is given by two landmarks in the space, whose positions are known.
 * The robot can measure the direct distance to both the landmarks, for instance
 * through visual localization techniques.
 *
 * @param T Numeric scalar type
 * @param CovarianceBase Class template to determine the covariance representation
 *                       (as covariance matrix (StandardBase) or as lower-triangular
 *                       coveriace square root (SquareRootBase))
 */
template<typename T, template<class> class CovarianceBase = Kalman::StandardBase>
class PositionMeasurementModel : public Kalman::LinearizedMeasurementModel<State<T>, PositionMeasurement<T>, CovarianceBase>
{
public:
    //! State type shortcut definition
    typedef  KalmanExamples::Robot1::State<T> S;
    
    //! Measurement type shortcut definition
    typedef  KalmanExamples::Robot1::PositionMeasurement<T> M;
    
    /**
     * @brief Constructor
     *
     * @param landmark1x The x-position of landmark 1
     * @param landmark1y The y-position of landmark 1
     * @param landmark2x The x-position of landmark 2
     * @param landmark2y The y-position of landmark 2
     */
    PositionMeasurementModel(T landmark1x, T landmark1y, T landmark2x, T landmark2y)
    {
        // Save landmark positions
        landmark1 << landmark1x, landmark1y;
        landmark2 << landmark2x, landmark2y;
        
        // Setup noise jacobian. As this one is static, we can define it once
        // and do not need to update it dynamically
        this->V.setIdentity();
    }
    
    /**
     * @brief Definition of (possibly non-linear) measurement function
     *
     * This function maps the system state to the measurement that is expected
     * to be received from the sensor assuming the system is currently in the
     * estimated state.
     *
     * @param [in] x The system state in current time-step
     * @returns The (predicted) sensor measurement for the system state
     */
    M h(const S& x) const
    {
        M measurement;
        
        // Robot position as (x,y)-vector
        // This uses the Eigen template method to get the first 2 elements of the vector
        Kalman::Vector<T, 2> position = x.template head<2>();
        
        // Distance of robot to landmark 1
        Kalman::Vector<T, 2> delta1 = position - landmark1;
        measurement.d1() = std::sqrt( delta1.dot(delta1) );
        
        // Distance of robot to landmark 2
        Kalman::Vector<T, 2> delta2 = position - landmark2;
        measurement.d2() = std::sqrt( delta2.dot(delta2) );
        
        return measurement;
    }
    
protected:
    //! Position of landmark 1 given as (x,y)-measurement
    Kalman::Vector<T, 2> landmark1;
    
    //! Position of landmark 2 given as (x,y)-measurement
    Kalman::Vector<T, 2> landmark2;

protected:
    
    /**
     * @brief Update jacobian matrices for the system state transition function using current state
     *
     * This will re-compute the (state-dependent) elements of the jacobian matrices
     * to linearize the non-linear measurement function \f$h(x)\f$ around the
     * current state \f$x\f$.
     *
     * @note This is only needed when implementing a LinearizedSystemModel,
     *       for usage with an ExtendedKalmanFilter or SquareRootExtendedKalmanFilter.
     *       When using a fully non-linear filter such as the UnscentedKalmanFilter
     *       or its square-root form then this is not needed.
     *
     * @param x The current system state around which to linearize
     * @param u The current system control input
     */
    void updateJacobians( const S& x )
    {
        // H = dh/dx (Jacobian of measurement function w.r.t. the state)
        this->H.setZero();
        
        // Robot position as (x,y)-vector
        // This uses the Eigen template method to get the first 2 elements of the vector
        Kalman::Vector<T, 2> position = x.template head<2>();
        
        // Distance of robot to landmark 1
        Kalman::Vector<T, 2> delta1 = position - landmark1;
        
        // Distance of robot to landmark 2
        Kalman::Vector<T, 2> delta2 = position - landmark2;
        
        // Distances
        T d1 = std::sqrt( delta1.dot(delta1) );
        T d2 = std::sqrt( delta2.dot(delta2) );
        
        // partial derivative of meas.d1() w.r.t. x.x()
        this->H( M::D1, S::X ) = delta1[0] / d1;
        // partial derivative of meas.d1() w.r.t. x.y()
        this->H( M::D1, S::Y ) = delta1[1] / d1;
        
        // partial derivative of meas.d1() w.r.t. x.x()
        this->H( M::D2, S::X ) = delta2[0] / d2;
        // partial derivative of meas.d1() w.r.t. x.y()
        this->H( M::D2, S::Y ) = delta2[1] / d2;
    }
};

} // namespace Robot
} // namespace KalmanExamples

#endif