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

#include "../header/gate.h"

using namespace arma;

/*********************************************************/
Gate::Gate(mat_ptr mat, set_mat matbw) : mat{mat}, matbw{matbw} {}

/*********************************************************/
Gate::Gate(str path) {
  mat = std::make_shared<sp_cx_mat>();
  mat->load(path, arma_binary);
}

/*********************************************************/
void Gate::save(str path){
  mat->save(path, arma_binary);
}

/*********************************************************/
str Gate::__str__() {
  sstr out;
  size_t size = log2(mat->n_rows);
  if (size == 1)
    out <<"1 qubit gate" << std::endl;
  else 
    out << size <<" qubits gate" << std::endl;

  for (auto i = mat->begin(); i != mat->end(); ++i) {
    auto aux = utility::cx_to_str(*i);
    out << "(" << i.row() << ", " << i.col() << std::left
        << std::setw(10) << ")"
        << (aux == ""? "1" : aux)  << std::endl;
  }
  return out.str();
}

/*********************************************************/
mat_ptr& Gate::get_mat() {
  return mat;
}

set& Gate::get_matbw(size_t i) {
  return matbw[i];
}

/*********************************************************/
Gate Gate::make_gate(vec_complex matrix) {
  if (matrix.size() != 4) {
    sstr err;
    err << "\'matrix\' argument must have exactly 4 elements: "
        << "[u00, u01, u10, u11]";
    throw std::invalid_argument{err.str()};
  }

  set_mat matbw;
  matbw[0].insert(std::make_pair(matrix[0], 0));
  matbw[0].insert(std::make_pair(matrix[2], 1));
  matbw[1].insert(std::make_pair(matrix[1], 0));
  matbw[1].insert(std::make_pair(matrix[3], 1));
  return Gate{std::make_shared<sp_cx_mat>(cx_mat{{{matrix[0], matrix[1]},
                                             {matrix[2], matrix[3]}}}),
         matbw};
}

/*********************************************************/
Gate Gate::make_mgate(size_t size, 
                  vec_size_t row,
                  vec_size_t col,
                 vec_complex value) {
  if (row.size() != col.size()
      or row.size() != value.size()
      or col.size() != value.size()) {
    sstr err;
    err << "Arguments \'row\', \'col\' and \'value\' must have the same size";
    throw std::invalid_argument{err.str()};
  }

  auto sizem = 1ul << size;
  auto mp = std::make_shared<sp_cx_mat>(sizem, sizem);
  auto &m = *mp;

  set_mat matbw;

  for (size_t i = 0; i < row.size(); i++) {
    m(row[i], col[i]) = value[i];
    matbw[col[i]].insert(std::make_pair(value[i], row[i]));
  }

  return Gate{mp, matbw};
}

Gate Gate::make_cgate(str gates,
                 vec_size_t control) {
  if (control.size() == 0) {
    sstr err;
    err << "\'control\' argument must have at least one item";
    throw std::invalid_argument{err.str()};
  }

  size_t size = gates.size();

  for (auto& i : control) {
    if (i >= size) {
      sstr err;
      err << "Items in \'control\' should be in the range of 0 to "
          << (size-1);
      throw std::invalid_argument{err.str()};
    }
  }

  size_t x = 0;
  size_t z = 0;
  for (size_t i = 0; i < size; i++) {
    if (gates[i] == 'X') {
      x |= 1ul << (size-i-1);
    } else if (gates[i] == 'Z') {
      z |= 1ul << (size-i-1);
    } else if (gates[i] == 'I') {
      continue;
    } else {
      sstr err;
      err << "Argument \'gates\' must have only \'X\', \'Z\' and \'I\'";
      throw std::invalid_argument{err.str()};
    }
  }

  auto parity = [&](size_t x) {
    for (int i = 32; i > 0; i /= 2)
      x ^= x >> i;
    return x & 1; 
  };

  auto cmp = std::make_shared<sp_cx_mat>(1ul << size, 1ul << size);
  auto& cm = *cmp;
  
  set_mat matbw;

  for (size_t i = 0; i < (1ul << size); i++) {
    bool cond = true;
    for (size_t k = 0; k < control.size(); k++)
      cond = cond and ((1ul << (size-control[k]-1)) & i);

    if (cond) {
      size_t row = (i ^ x);
      cm(row, i) = pow(-1, parity(i & z));
      matbw[i].insert(std::make_pair(pow(-1, parity(i & z)), row));
    } else {
      cm(i,i) = 1;
      matbw[i].insert(std::make_pair(1, i));
    }
  }

  return Gate{cmp, matbw};
}

Gate Gate::make_fgate(PyObject* func,
                   size_t size,
                PyObject* iterator) {

  auto mp = std::make_shared<sp_cx_mat>(1ul << size, 1ul << size);
  auto &m = *mp;

  if (iterator == Py_None) {
    PyObject *builtins = PyEval_GetBuiltins(); 
    PyObject *range = PyDict_GetItemString(builtins , "range");
    iterator = PyEval_CallFunction(range, "i", 1ul << size);
  }

  auto* it = PyObject_GetIter(iterator);

  set_mat matbw;
  
  PyObject* pyj;
  while ((pyj = PyIter_Next(it))) {
    auto* arg = PyTuple_Pack(1, pyj);
    auto* pyi = PyObject_CallObject(func,   arg);
    
    auto i = PyLong_AsSize_t(pyi);
    auto j = PyLong_AsSize_t(pyj);

    m(i, j) = 1;
    matbw[j].insert(std::make_pair(1, i));
   
    Py_DECREF(pyi);
    Py_DECREF(pyj);
  }

  Py_DECREF(it);

  return Gate{mp, matbw};
}

