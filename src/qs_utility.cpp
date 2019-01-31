/* MIT License
 * 
 * Copyright (c) 2019 Evandro Chagas Ribeiro da Rosa
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

#include "../header/qsystem.h"
#include <eigen3/unsupported/Eigen/KroneckerProduct>
#include <iomanip>
#include <sstream>

QSystem::QSystem(size_t nqbits, size_t seed, Gate& gate, std::string state) :
  gate{gate}, size{nqbits}, state{state}, ops{new char[size]}, syncc{true},
  qbits{1l << size, state == "mix" ? 1l << size : 1},
  bits{new Bit[size]()}, an_size{0}, an_ops{nullptr}, an_bits{nullptr}
{
  if (state != "mix" and state != "pure") 
    throw std::invalid_argument{"Argument \'state\' must be \"pure\" or \"mix\", not \""
      + state + "\"."};
  qbits.coeffRef(0,0) = 1;
  std::memset(ops, 'I', size*sizeof(char));
  std::srand(seed);
}

QSystem::~QSystem() {
  delete[] ops;
  delete[] bits;
  if (an_ops) delete[] an_ops;
  if (an_bits) delete[] an_bits;
}

std::string QSystem::__str__() {
  auto to_bits = [&](size_t i) {
    std::string sbits{'|'};
    for (size_t j = 0; j < size; j++)
      sbits += i & 1ul << (size+an_size-j-1)? '1' : '0';
    sbits += an_size == 0? ">" : ">|";
    for (size_t j = size; j < size+an_size; j++)
      sbits += i & 1ul << (size+an_size-j-1)? '1' : '0';
    sbits += an_size == 0? "" : ">";
    return sbits;    
  };

  auto cx_to_str = [&](cx i) {
    std::stringstream ss;
    if (fabs(i.imag()) < 1e-14) {
      ss << std::showpos << std::fixed
         << std::setprecision(3)  << i.real() << "       ";
    } else if (fabs(i.real()) < 1e-14) {
      ss << std::showpos << std::fixed
         << std::setprecision(3) << std::setw(12) << i.imag() << 'i';
    } else {
      ss << std::showpos << std::fixed << std::setprecision(3) << i.real()
         << i.imag() << 'i';
    }
    return ss.str();
  };

  if (not syncc) sync();
  std::stringstream out;
  if (state == "pure") {
    for (it_mat i(qbits, 0); i; ++i) {
      if (abs(i.value()) < 1e-14) continue; 
      out << cx_to_str(i.value()) << to_bits(i.row()) << '\n';
    }
  } else if (state == "mix") {
    for (auto k = 0l; k < qbits.outerSize(); ++k) {
      for (it_mat i(qbits, k); i; ++i) {
        auto aux = cx_to_str(i.value());
        out << "(" << i.row() << ", " << i.col() << ")    " <<
          (aux == ""? "1" : aux)  << std::endl;
      }
    }
  }
  return out.str();
}

size_t QSystem::get_size() {
  return size;
}

std::vector<int> QSystem::get_bits() {
  std::vector<int> vec;
  for (size_t i = 0; i < size; i++)
    vec.push_back(bits[i]);
  return vec;
}

size_t QSystem::get_an_size() {
  return an_size;
}

std::vector<int> QSystem::get_an_bits() {
  std::vector<int> vec;
  for (size_t i = 0; i < an_size; i++)
    vec.push_back(an_bits[i]);
  return vec;
}

PyObject* QSystem::get_qbits() {
  if (not syncc) sync();
  PyObject* top_tuple = PyTuple_New(2);
  PyObject* rc_tuple = PyTuple_New(2);
  PyObject* val = PyList_New(qbits.nonZeros());
  PyObject* row = PyList_New(qbits.nonZeros());
  PyObject* col = PyList_New(qbits.nonZeros());
  size_t v = 0;
  size_t r = 0;
  size_t c = 0;
  for (auto k = 0l; k < qbits.outerSize(); ++k) {
    for (it_mat i(qbits, k); i; ++i) {
      PyList_SetItem(val, v++, PyComplex_FromDoubles(i.value().real(),
                                                     i.value().imag()));
      PyList_SetItem(row, r++, PyLong_FromLong(i.row()));
      PyList_SetItem(col, c++, PyLong_FromLong(i.col()));
    } 
  }

  PyTuple_SetItem(rc_tuple, 0, row);
  PyTuple_SetItem(rc_tuple, 1, col);
  PyTuple_SetItem(top_tuple, 0, val);
  PyTuple_SetItem(top_tuple, 1, rc_tuple);

  PyObject* size_tuple = PyTuple_New(2);
  PyTuple_SetItem(size_tuple, 0, PyLong_FromLong(1l << (size+an_size)));
  PyTuple_SetItem(size_tuple, 1, PyLong_FromLong(state == "pure"? 1l :
                                                 1l << (size+an_size)));

  PyObject* result = PyTuple_New(2);
  PyTuple_SetItem(result, 0, top_tuple);
  PyTuple_SetItem(result, 1, size_tuple);

  return result;
}

void QSystem::change_to(std::string state) {
  if (state != "mix" and state != "pure") 
    throw std::invalid_argument{"Argument \'state\' must be \"pure\" or \"mix\", not \""
      + state + "\"."};

  if (state == this->state) 
    return;

  if (state == "mix") {
    qbits = qbits*qbits.adjoint();
  } else if (state == "pure") {
    sp_cx_mat nqbits{1l << (size+an_size), 1};
    sp_cx_vec diag = qbits.diagonal().sparseView();
    for (it_vec i(diag); i; ++i)
      nqbits.coeffRef(i.index(), 0) = sqrt(i.value().real());
    qbits = nqbits;
  }
  
  this->state = state;
}

std::string QSystem::get_state() {
  return state;
}

sp_cx_mat QSystem::make_gate(sp_cx_mat gate, size_t qbit) {
  sp_cx_mat m;
  if (qbit == 0) {
    size_t eyesize = 1l << (size+an_size-1);
    m = kron(gate, eye(eyesize, eyesize));
  } else if (qbit == size+an_size-1) {
    size_t eyesize = 1l << (size+an_size-1);
    m = kron(eye(eyesize, eyesize), gate);
  } else {
    size_t eyesize = 1l << qbit;
    m = kron(eye(eyesize, eyesize), gate);
    eyesize = 1l << (size+an_size-qbit-1);
    m = kron(m, eye(eyesize, eyesize));
  }
  return m;
}

sp_cx_mat QSystem::kron(const sp_cx_mat& a, const sp_cx_mat& b) {
  return kroneckerProduct(a,b).eval();
}

sp_cx_mat QSystem::eye(long int rows, long int cols) {
  sp_cx_mat m(rows, cols);
  m.setIdentity();
  return m;
}

