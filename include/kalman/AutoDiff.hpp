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
#ifndef KALMAN_AUTODIFF_HPP
#define KALMAN_AUTODIFF_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/AutoDiff>
#include "Matrix.hpp"

namespace Kalman {

    template<typename T, template<typename> class Input, template<typename> class Output>
    class AutoDiffJacobian {
        static_assert(Input<T>::ColsAtCompileTime == 1, "Input must be column vector");
        static_assert(Output<T>::ColsAtCompileTime == 1, "Output must be column vector");

    public:
        typedef Eigen::Matrix<T, Output<T>::RowsAtCompileTime, Input<T>::RowsAtCompileTime> Jacobian;

        typedef Eigen::Matrix<T, Input<T>::RowsAtCompileTime, 1> Gradient;
        typedef Eigen::AutoDiffScalar<Gradient> ActiveScalar;

        typedef Vector<ActiveScalar, Input<T>::RowsAtCompileTime> ActiveInput;
        typedef Vector<ActiveScalar, Output<T>::RowsAtCompileTime> ActiveOutput;

        typedef std::function<void(const ActiveInput&, ActiveOutput&)> Function;

    public:

        AutoDiffJacobian(const Function& func) : function(func) {

        }

        void operator()(const Input<T>& input, Output<T>& output, Jacobian& jacobian)
        {
            ActiveInput in = input.template cast<ActiveScalar>();
            ActiveOutput out;
            for(int i = 0; i < input.rows(); i++)
            {
                in[i].derivatives() = Gradient::Unit(input.rows(), i);
            }

            if(Output<T>::RowsAtCompileTime == Dynamic) {
                out.resize(output.rows());
                for(size_t i = 0; i < input.rows(); i++)
                {
                    out[i].derivatives().resize(input.rows());
                }
            }

            function(in, out);

            if(Input<T>::RowsAtCompileTime == Dynamic || Output<T>::RowsAtCompileTime == Dynamic)
            {
                jacobian.resize(out.rows(), in.rows());
            }

            for(int i = 0; i < out.rows(); i++)
            {
                output[i]       = out[i].value();
                jacobian.row(i) = out[i].derivatives();
            }
        }
    protected:
        const Function& function;
    };

}

#endif //KALMAN_AUTODIFF_HPP
