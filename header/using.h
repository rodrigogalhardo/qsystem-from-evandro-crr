/* MIT License
 * 
 * Copyright (c) 2019 Evandro Chagas Ribeiro da Rosa <ev.crr97@gmail.com>
 * Copyright (c) 2019 Bruno Gouvêa Taketani <b.taketani@ufsc.br>
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
#include <cstring>
#include <complex>
#include <utility>
#include <sstream>
#include <fstream>
#include <memory>
#include <Python.h>
#include <armadillo>
#include <iomanip>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/serialization/boost_unordered_map.hpp>
#include <boost/serialization/boost_unordered_set.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

using size_t = long unsigned;
using complex = std::complex<double>;
using str = std::string;
using vec_complex = std::vector<std::complex<double>>;
using vec_size_t = std::vector<size_t>;
using vec_str = std::vector<str>;
using vec_int = std::vector<int>;
using vec_float = std::vector<double>;
using cnot_pair = std::pair<size_t, vec_size_t>;
using cph_tuple = std::tuple<complex, size_t, vec_size_t>;
using u3_tuple = std::tuple<double, double, double>;
using cut_pair = std::pair<size_t, size_t>;
using r_pair = std::pair<char, double>;
using sstr = std::stringstream;
using mat_ptr = std::shared_ptr<arma::sp_cx_mat>;
using py_obj = PyObject*;
using py_function = PyObject*;
using py_iterator = PyObject*;
using dict = boost::unordered_map<size_t, complex>;
using set = boost::unordered_set<std::pair<complex, size_t>>;
using set_mat = boost::unordered_map<size_t, set>;

namespace utility {
  str cx_to_str(complex i, bool use_sqrt = true);
  str to_bits(size_t i, size_t qsize, size_t asize);
}
