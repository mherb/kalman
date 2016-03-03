// The MIT License (MIT)
//
// Copyright (c) 2015 Markus Herb
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
#ifndef KALMAN_MATRIX_HPP_
#define KALMAN_MATRIX_HPP_

#include <cmath>

#include <Eigen/Dense>

#define KALMAN_VECTOR(NAME, T, N)                                                       \
    typedef Kalman::Vector<T, N> Base;                                                  \
    using typename Base::Scalar;                                                        \
    using Base::RowsAtCompileTime;                                                      \
    using Base::ColsAtCompileTime;                                                      \
    using Base::SizeAtCompileTime;                                                      \
                                                                                        \
    NAME(void) : Kalman::Vector<T, N>() {}                                              \
                                                                                        \
    template<typename OtherDerived>                                                     \
    NAME(const Eigen::MatrixBase<OtherDerived>& other) : Kalman::Vector<T, N>(other) {} \
                                                                                        \
    template<typename OtherDerived>                                                     \
    NAME& operator= (const Eigen::MatrixBase <OtherDerived>& other)                     \
    {                                                                                   \
        this->Base::operator=(other);                                                   \
        return *this;                                                                   \
    }

namespace Kalman {
    const int Dynamic = Eigen::Dynamic;

    /**
     * @class Kalman::Matrix
     * @brief Template type for matrices
     * @param T The numeric scalar type
     * @param rows The number of rows
     * @param cols The number of columns
     */
    template<typename T, int rows, int cols>
    using Matrix = Eigen::Matrix<T, rows, cols>;
    
    /**
     * @brief Template type for vectors
     * @param T The numeric scalar type
     * @param N The vector dimension
     */
    template<typename T, int N>
    class Vector : public Matrix<T, N, 1>
    {
    public:
        //! Matrix base type
        typedef Matrix<T, N, 1> Base;

        using typename Base::Scalar;
        using Base::RowsAtCompileTime;
        using Base::ColsAtCompileTime;
        using Base::SizeAtCompileTime;

        Vector(void) : Matrix<T, N, 1>() {}
        
        /**
         * @brief Copy constructor
         */
        template<typename OtherDerived>
        Vector(const Eigen::MatrixBase<OtherDerived>& other)
            : Matrix<T, N, 1>(other)
        { }
        /**
         * @brief Copy assignment constructor
         */
        template<typename OtherDerived>
        Vector& operator= (const Eigen::MatrixBase <OtherDerived>& other)
        {
            this->Base::operator=(other);
            return *this;
        }
    };
    
    /**
     * @brief Cholesky square root decomposition of a symmetric positive-definite matrix
     * @param _MatrixType The matrix type
     * @param _UpLo Square root form (Eigen::Lower or Eigen::Upper)
     */
    template<typename _MatrixType, int _UpLo = Eigen::Lower>
    class Cholesky : public Eigen::LLT< _MatrixType, _UpLo >
    {
    public:
        Cholesky() : Eigen::LLT< _MatrixType, _UpLo >() {}
        
        /**
         * @brief Construct cholesky square root decomposition from matrix
         * @param m The matrix to be decomposed
         */
        Cholesky(const _MatrixType& m ) : Eigen::LLT< _MatrixType, _UpLo >(m) {}
        
        /**
         * @brief Set decomposition to identity
         */
        Cholesky& setIdentity()
        {
            this->m_matrix.setIdentity();
            this->m_isInitialized = true;
            return *this;
        }
        
        /**
         * @brief Check whether the decomposed matrix is the identity matrix
         */
        bool isIdentity() const
        {
            eigen_assert(this->m_isInitialized && "LLT is not initialized.");
            return this->m_matrix.isIdentity();
        }
        
        /**
         * @brief Set lower triangular part of the decomposition
         * @param matrix The lower part stored in a full matrix
         */
        template<typename Derived>
        Cholesky& setL(const Eigen::MatrixBase <Derived>& matrix)
        {
            this->m_matrix = matrix.template triangularView<Eigen::Lower>();
            this->m_isInitialized = true;
            return *this;
        }
        
        /**
         * @brief Set upper triangular part of the decomposition
         * @param matrix The upper part stored in a full matrix
         */
        template<typename Derived>
        Cholesky& setU(const Eigen::MatrixBase <Derived>& matrix)
        {
            this->m_matrix = matrix.template triangularView<Eigen::Upper>().adjoint();
            this->m_isInitialized = true;
            return *this;
        }
    };
}

#endif
