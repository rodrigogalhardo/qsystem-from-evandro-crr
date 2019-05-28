/* MIT License
 * 
 * Copyright (c) 2019 Bruno Gouvêa Taketani <b.taketani@ufsc.br>
 * Copyright (c) 2019 Evandro Chagas Ribeiro da Rosa <ev.crr97@gmail.com>
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */     

#pragma once
#include <vector>
#include <string>
#include <complex>
#include <utility>
#include <sstream>
#include <memory>
#include <Python.h>
#include <armadillo>
using size_t = long unsigned;
using complex = std::complex<double>;
using vec_complex = std::vector<std::complex<double>>;
using vec_size_t = std::vector<size_t>;
using vec_str = std::vector<std::string>;
using vec_int = std::vector<int>;
using vec_float = std::vector<double>;
using cnot_pair = std::pair<size_t, vec_size_t>;
using cph_tuple = std::tuple<complex, size_t, vec_size_t>;
using u3_tuple = std::tuple<double, double, double>;
using cut_pair = std::pair<size_t, size_t>;
using r_pair = std::pair<char, double>;
using sstr = std::stringstream;
using str = std::string;
using mat_ptr = std::shared_ptr<arma::sp_cx_mat>;
using py_obj = PyObject*;
using py_function = PyObject*;
using py_iterator = PyObject*;
