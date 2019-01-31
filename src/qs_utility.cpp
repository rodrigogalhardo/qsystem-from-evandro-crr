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

/*
QSystem::QSystem(std::string path, size_t seed, Gate& gate) :
  gate{gate}, syncc{true}, an_size{0}, an_ops{nullptr}, an_bits{nullptr}
{
  qbits.load(path, arma_binary);
  size = log2(qbits.n_rows);
  state = qbits.n_cols > 1 ? "mix" : "pure";
  ops = new char[size];
  bits = new Bit[size]();
  std::memset(ops, 'I', size*sizeof(char));
  std::srand(seed);
}
*/

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

  auto cx_to_str = [&](std::complex<double> i) {
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

/*
PyObject* QSystem::get_qbits() {
  if (not syncc) sync();
  qbits.sync();

  PyObject* csc_tuple = PyTuple_New(3);
  PyObject* val = PyList_New(qbits.n_nonzero);
  PyObject* row_ind = PyList_New(qbits.n_nonzero);
  for (size_t i = 0; i < qbits.n_nonzero; i++) {
    PyList_SetItem(val, i, PyComplex_FromDoubles(qbits.values[i].real(),
                                                 qbits.values[i].imag()));
    PyList_SetItem(row_ind, i, PyLong_FromLong(qbits.row_indices[i]));
  }
  PyTuple_SetItem(csc_tuple, 0, val);
  PyTuple_SetItem(csc_tuple, 1, row_ind);

  PyObject* col_ptr = PyList_New(qbits.n_cols+1);
  for (size_t i = 0; i < qbits.n_cols+1; i++) 
    PyList_SetItem(col_ptr, i, PyLong_FromLong(qbits.col_ptrs[i]));
  PyTuple_SetItem(csc_tuple, 2, col_ptr);

  PyObject* size_tuple = PyTuple_New(2);
  PyTuple_SetItem(size_tuple, 0, PyLong_FromLong(qbits.n_rows));
  PyTuple_SetItem(size_tuple, 1, PyLong_FromLong(qbits.n_cols));

  PyObject* result = PyTuple_New(2);

  PyTuple_SetItem(result, 0, csc_tuple);
  PyTuple_SetItem(result, 1, size_tuple);

  return result;
}

void QSystem::set_qbits(vec_size row_ind,
                        vec_size col_ptr,
                          vec_cx values,
                          size_t nqbits,
                     std::string state) {
  if (not syncc) clar();

  qbits = sp_cx_mat(conv_to<uvec>::from(row_ind),
                    conv_to<uvec>::from(col_ptr),
                    cx_vec(values),
                    1ul << nqbits,
                    state == "pure"? 1ul : 1ul << nqbits);
                    
  this->state = state;
  size = nqbits;
}
*/

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

/*
void QSystem::save(std::string path) {
  if (not syncc) sync();
  qbits.save(path, arma_binary);
}
*/

