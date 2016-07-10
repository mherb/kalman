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

#ifndef KALMAN_DEPRECATED_HPP_
#define KALMAN_DEPRECATED_HPP_

// Adapted from http://stackoverflow.com/questions/295120/c-mark-as-deprecated/21265197#21265197


#ifndef DEPRECATED
    #ifdef __GNUC__
        #define DEPRECATED __attribute__((deprecated))
    #elif defined(_MSC_VER)
        #define DEPRECATED __declspec(deprecated)
    #else
        #pragma message("WARNING: You need to implement DEPRECATED for this compiler")
        #define DEPRECATED
    #endif
#endif

#endif //KALMAN_DEPRECATED_HPP_
