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
#include <algorithm>

using namespace arma;

/******************************************************/
void QSystem::evol(std::string gate,
                        size_t qbit, 
                        size_t count,
                          bool inver) {
  if (qbit >= size()) {
      sstr err;
      err << "\'qbit\' argument should be in the range of 0 to "
          << (size()-1);
      throw std::invalid_argument{err.str()};
  } else if (count == 0 and qbit+count <= size()) {
      sstr err;
      err << "\'cout\' argument should be greater than 0 "
          << "and \'qbit+count\' suld be in the range of 0 to "
          << size();
      throw std::invalid_argument{err.str()};
  }

  if (gate.size() > 1) {
    auto size_n = log2(gates.mget(gate).n_rows);
    sync(qbit, qbit+count*size_n);
    for (size_t i = 0; i < count; i++) {
      size_t index = qbit+i*gate.size();
      fill(Gate_aux::GATE_N, index, size_n);
      ops(index).data = gate;
      ops(index).inver = inver;

    }
  } else {
    sync(qbit, qbit+count);
    for (size_t i = 0; i < count; i++) {
      ops(qbit+i).tag = Gate_aux::GATE_1;
      ops(qbit+i).data = gate[0];
      ops(qbit+i).inver = inver;
    }
  }

  _sync = false;
}

/******************************************************/
void QSystem::cnot(size_t target, vec_size_t control) {
  if (target >= size()) {
      sstr err;
      err << "\'target\' argument should be in the range of 0 to "
          << (size()-1);
      throw std::invalid_argument{err.str()};
  } else if (control.size() == 0) {
    sstr err;
    err << "\'control\' argument must have at least one item";
    throw std::invalid_argument{err.str()};
  }
  for (auto& i : control) {
    if (i >= size()) {
      sstr err;
      err << "Items in \'control\' should be in the range of 0 to "
          << (size()-1);
      throw std::invalid_argument{err.str()};
    }
  }

  auto [size_n, minq] = cut(target, control);
  fill(Gate_aux::CNOT, minq, size_n);
  ops(minq).data = cnot_pair{target, control};
}

/******************************************************/
void QSystem::cphase(complex phase, size_t target, vec_size_t control) {
  if (target >= size()) {
      sstr err;
      err << "\'target\' argument should be in the range of 0 to "
          << (size()-1);
      throw std::invalid_argument{err.str()};
  } else if (std::abs(std::abs(phase) - 1.0) > 1e-14) {
    sstr err;
    err << "abs(phase) must be equal to 1";
    throw std::invalid_argument{err.str()};
  } else if (control.size() == 0) {
    sstr err;
    err << "\'control\' argument must have at least one item";
    throw std::invalid_argument{err.str()};
  }
  for (auto& i : control) {
    if (i >= size()) {
      sstr err;
      err << "Items in \'control\' should be in the range of 0 to "
          << (size()-1);
      throw std::invalid_argument{err.str()};
    }
  }

  auto [size_n, minq] = cut(target, control);
  fill(Gate_aux::CPHASE, minq, size_n);
  ops(minq).data = cph_tuple{phase, target, control};
}

/******************************************************/
void QSystem::swap(size_t qbit_a, size_t qbit_b) {
  if (qbit_a >= size() or qbit_b >= size()) {
      sstr err;
      err << "Arguments \'qbit_a\' and \'qbit_b\' should be in the "
          << "range of 0 to " << (size()-1);
      throw std::invalid_argument{err.str()};
  }

  if (qbit_a == qbit_b) return;
  size_t a = qbit_a < qbit_b? qbit_a :  qbit_b;
  size_t b = qbit_a > qbit_b? qbit_a :  qbit_b;
  fill(Gate_aux::SWAP, a, b-a+1);
}

/******************************************************/
void QSystem::qft(size_t qbegin, size_t qend, bool inver) {
  if (qbegin >= size() or qend > size() or qbegin <= qend) {
      sstr err;
      err << "\'qbegin\' argument should be in the "
          << "range of 0 to " << (size()-1)
          << " and argument \'qend\' should be greater than 0 "
          << "and in the range of 0 to " << size();
      throw std::invalid_argument{err.str()};
  }

  fill(Gate_aux::QFT, qbegin, qend-qbegin);
  ops(qbegin).inver = inver;
}

/******************************************************/
void QSystem::sync() {
  
  if (_sync) return;

  sp_cx_mat evolm;
  
  evolm = get_gate(ops(0));

  for (size_t i = ops(0).size; i < size(); i += ops(i).size) {
    evolm = kron(evolm, get_gate(ops(i)));
  }

  if (_state == "vector")
    qbits = evolm*qbits;
  else if (_state == "matrix")
    qbits = evolm*qbits*evolm.t();

  delete[] _ops;
  _ops = new Gate_aux[_size]();
  if (an_ops) {
    delete[] an_ops;
    an_ops = new Gate_aux[an_size]();
  }

  _sync = true;
}

/******************************************************/
void QSystem::sync(size_t qbegin, size_t qend) {
  for (size_t i = qbegin; i < qend; i++) {
    if (ops(i).busy()) {
      sync();
      break;
    } 
  }
}

/******************************************************/
void QSystem::clar() {
  delete _ops;
  _ops = new Gate_aux[_size]();
  if (an_ops) {
    delete an_ops;
    an_ops = new Gate_aux[an_size]();
  }

  _sync = true;
}

