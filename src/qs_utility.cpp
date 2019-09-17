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

#include "../header/qsystem.h"

using namespace arma;

/******************************************************/
QSystem::QSystem(size_t num_qbits,
                  size_t seed,
             std::string representation,
                  size_t init) :
  _size{num_qbits},
  _representation{representation},
  _sync{true},
  _bits{new Bit[num_qbits]()}, 
  an_size{0},
  an_ops{nullptr},
  an_bits{nullptr}
{
  valid_state_str(representation);
  valid_init(init, num_qbits);
  if (representation == "bitwise") {
    bwqbits[init] = 1;
    _ops = nullptr;
  } else {
    qbits = sp_cx_mat{1lu << num_qbits, representation == "matrix" ? 1lu << num_qbits : 1};
    qbits(init, representation == "matrix"? init : 0) = 1;
    _ops = new Gate_aux[num_qbits]();
  }

  std::srand(seed);
}


/******************************************************/
QSystem::~QSystem() {
  if (_ops) delete[] _ops;
  delete[] _bits;
  if (an_ops) delete[] an_ops;
  if (an_bits) delete[] an_bits;
}

/******************************************************/
QSystem::Gate_aux::Gate_aux() : tag{GATE_1}, data{'I'}, size{1}, invert{false} {}

/******************************************************/
QSystem::Gate_aux::~Gate_aux() {}

/******************************************************/
bool QSystem::Gate_aux::busy() {
  return not(tag == GATE_1 and std::get<char>(data) == 'I');
}

/******************************************************/
std::string QSystem::__str__() {
  sync();
  sstr out;
  if (representation() == "vector") {
    for (auto i = qbits.begin(); i != qbits.end(); ++i) {
      if (abs((cx_double)*i) < 1e-14) continue; 
      out << utility::cx_to_str(*i)
          << utility::to_bits(i.row(), _size, an_size) << std::endl;
    }
  } else if (representation() == "matrix") {
    for (auto i = qbits.begin(); i != qbits.end(); ++i) {
      auto aux = utility::cx_to_str(*i, false);
      out << "(" << i.row() << ", " << i.col() << std::left
          << std::setw(10) << ")"
          << (aux == ""? "1" : aux)  << std::endl;
    }
  } else if (representation() == "bitwise") {
    for (auto &ket : bwqbits) 
      out << utility::cx_to_str(ket.second)
          << utility::to_bits(ket.first, _size, an_size) 
          << std::endl;
  }
  return out.str();
}

/******************************************************/
size_t QSystem::size() {
  return _size+an_size;
}

/******************************************************/
PyObject* QSystem::get_qbits() {
  valid_not_bw();
  sync();
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

/******************************************************/
void QSystem::set_qbits(vec_size_t row_ind,
                        vec_size_t col_ptr,
                       vec_complex values,
                            size_t num_qbits,
                       std::string representation) {
  qbits = sp_cx_mat(conv_to<uvec>::from(row_ind),
                    conv_to<uvec>::from(col_ptr),
                    cx_vec(values),
                    1ul << num_qbits,
                    representation == "vector"? 1ul : 1ul << num_qbits);
                    
  this->_representation = representation;
  _size = num_qbits;
  clear();
}

/******************************************************/
void QSystem::change_to(std::string new_state) {
  valid_state_str(new_state);

  if (new_state == representation()) 
    return;

  if (representation() == "vector") {
    if (new_state == "matrix") {
      qbits = qbits*qbits.t();
    } else if (new_state == "bitwise") {
      for (auto i = qbits.begin(); i != qbits.end(); ++i) {
        bwqbits[i.row()] = (complex) *i;
      } 
      qbits.zeros();
      delete[] _ops;
      if (an_ops) delete[] an_ops;
    }
  } else if (representation() == "matrix") {
    sstr err;
    err << "can not change the representation from \"matrix\"";
    throw std::runtime_error{err.str()};
  } else if (representation() == "bitwise") {
    qbits = sp_cx_mat{1ul << size(), 1};
    for (auto &ket : bwqbits) 
      qbits(ket.first, 1) = ket.second;
    
    if (new_state == "matrix") 
      qbits = qbits*qbits.t();

    _ops = new Gate_aux[_size]();
    an_ops = new Gate_aux[an_size]();
  }
  
  _representation = new_state;
}

/******************************************************/
std::string QSystem::representation() {
  return _representation;
}

/******************************************************/
void QSystem::save(std::string path) {

  std::ofstream outfile(path);
  if (representation() != "bitwise") {
    sync();
    outfile.write("QsMT", 4);
    qbits.save(outfile, arma_binary);
  } else {
    outfile.write("QsBW", 4);
    outfile << size();
    boost::archive::binary_oarchive oarchf(outfile);
    oarchf << bwqbits;
  }
  outfile.close();
}

void QSystem::load(std::string path) {
  std::ifstream infile(path);
  char tag[5];
  infile.read(tag, 4);
  tag[4] = '\0';

  if (std::strcmp("QsBW", tag) != 0) {
    qbits.load(infile, arma_binary);
    _size = log2(qbits.n_rows);
    _representation = qbits.n_cols > 1 ? "matrix" : "vector";
  } else {
    infile >> _size;
    boost::archive::binary_iarchive iarch(infile);
    iarch >> bwqbits;
    _representation = "bitwise";    
  }
  clear();
}

/******************************************************/
void QSystem::clear() {
  _sync = true;
  an_size = 0;
  if (an_ops) {
    delete[] an_ops;
    delete[] an_bits;
  }
  an_ops = nullptr; 
  an_bits = nullptr;

  delete[] _ops;
  delete[] _bits;
  if (representation() != "bitwise")
    _ops = new Gate_aux[_size]();
  _bits = new Bit[_size]();
}
