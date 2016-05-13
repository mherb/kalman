// The MIT License (MIT)
//
// Copyright (c) 2016 Markus Herb
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
#ifndef KALMAN_MODEL_HPP
#define KALMAN_MODEL_HPP

#include <type_traits>
#include "Types.hpp"

namespace Kalman
{
    /**
     * @brief Model Helper
     *
     * This is an internal helper class to resolve System and Measurement Model template function overloads
     */
    class Model
    {
    public:
        /**
         * @brief Model Covariance Getter (Template Overload)
         * @param Model The model type
         * @param T The numeric type of the covariance matrix
         * @param Type The state/measurement type
         * @param model The model instance
         * @return The covariance obtained from the model
         */
        template<class Model, typename T, class Type>
        static auto getCovariance(Model& model)
        -> decltype((void)model.template getCovariance<T>(), Covariance<T, Type>())
        {
            return model.template getCovariance<T>();
        }

        /**
         * @brief Model Covariance Getter (Cast Overload)
         * @param Model The model type
         * @param T The numeric type of the covariance matrix
         * @param Type The state/measurement type
         * @param model The model instance
         * @return The covariance obtained from the model
         */
        template<class Model, typename T, class Type>
        static auto getCovariance(Model& model)
        -> decltype((void)model.getCovariance(), Covariance<T, Type>())
        {
            return model.getCovariance().template cast<T>();
        }

        /**
         * @brief Model Covariance Square Root Getter (Template Overload)
         * @param Model The model type
         * @param T The numeric type of the covariance matrix
         * @param Type The state/measurement type
         * @param model The model instance
         * @return The covariance square root obtained from the model
         */
        template<class Model, typename T, class Type>
        static auto getCovarianceSquareRoot(Model& model)
        -> decltype((void)model.template getCovariancSquareRoote<T>(), CovarianceSquareRoot<T, Type>())
        {
            return model.template getCovarianceSquareRoot<T>();
        }

        /**
         * @brief Model Covariance Square Root Getter (Cast Overload)
         * @param Model The model type
         * @param T The numeric type of the covariance matrix
         * @param Type The state/measurement type
         * @param model The model instance
         * @return The covariance square root obtained from the model
         */
        template<class Model, typename T, class Type>
        static auto getCovarianceSquareRoot(Model& model)
        -> decltype((void)model.getCovarianceSquareRoot(), CovarianceSquareRoot<T, Type>())
        {
            return model.getCovarianceSquareRoot().template cast<T>();
        }

        /**
         * @brief Template Model Checking Trait
         * This trait may be used to check whether the model implements the new template-based
         * API for predict/measurement steps. This may be used for both system and measurement models
         * using different variadic arguments.
         *
         * @param Model The Model type
         * @param State The State type
         * @param Args... Additional Arguments
         */
        template<class Model, typename State, typename... Args>
        struct isTemplateModel {
            template<class _Model, typename _State, typename... _Args>
            static constexpr auto test(size_t)
            -> decltype(std::declval<_Model>().predict(std::declval<_State>(), std::declval<_Args>()..., std::declval<_State&>()), bool())
            {
                return true;
            }

            template<class _Model, typename _State, typename... _Args>
            static constexpr auto test(ssize_t)
            -> decltype(std::declval<_Model>().measure(std::declval<_State>(), std::declval<_Args&>()...), bool())
            {
                return true;
            }

            template<class _Model, typename _State, typename... _Args>
            static constexpr bool test(...)
            {
                return false;
            }

            static constexpr bool value = test<Model, State, Args...>(size_t());
        };

        /**
         * @brief Explicit Jacobian Trait
         * This trait may be used to check whether a model implements an explicit jacobian calculation
         * instead of relying on auto-diff for jacobian computation. This may be used for both system and measurement models
         * using different variadic arguments.
         *
         * @param Model The Model type
         * @param State The State type
         * @param Args... Additional Arguments
         */
        template<class Model, typename State, typename... Args>
        struct hasJacobian {
            template<class _Model, typename _State, typename... _Args>
            static constexpr auto test(size_t)
            -> decltype(std::declval<_Model>().getJacobian(std::declval<_State>(), std::declval<_Args>()...), bool())
            {
                return true;
            }

            template<class _Model, typename _State, typename... _Args>
            static constexpr bool test(...)
            {
                return false;
            }

            static constexpr bool value = test<Model, State, Args...>(size_t());
        };

        /**
         * @brief Legacy Model Checking Trait
         * This trait may be used to check whether the model implements the old inheritance-based
         * API for predict/measurement steps. This may be used for both system and measurement models
         * using different variadic arguments.
         *
         * @param Model The Model type
         * @param State The State type
         * @param Args... Additional Arguments
         */
        template<class Model, typename... Args>
        struct isLegacyModel {
            template<class _Model, typename... _Args>
            static constexpr auto test(size_t)
            -> decltype(std::declval<_Model>().f(std::declval<_Args>()...), bool())
            {
                return true;
            }

            template<class _Model, typename... _Args>
            static constexpr auto test(ssize_t)
            -> decltype(std::declval<_Model>().h(std::declval<_Args>()...), bool())
            {
                return true;
            }

            template<class _Model, typename _State, typename... _Args>
            static constexpr bool test(...)
            {
                return false;
            }

            static constexpr bool value = test<Model, Args...>(size_t());
        };
    };
}

#endif //KALMAN_MODEL_HPP
