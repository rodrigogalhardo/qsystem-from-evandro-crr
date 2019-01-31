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

void QSystem::evol(char gate, size_t qbit) {
  if (ops[qbit] != 'I') sync();

  ops[qbit] = gate;
  syncc = false;
}

void QSystem::evol(char gate, size_t qbegin, size_t qend) {
  for (size_t i = qbegin; i < qend; i++) 
    evol(gate, i);
}

void QSystem::evol(std::string gates) {
  if (gates.size() != size+an_size) 
    throw std::length_error{"String 'gates' must be " +
        std::to_string(size+an_size) + " characters long"};
  
  for (size_t i = 0; i < size; i++) 
    evol(gates[i], i);

  for (size_t i = 0; i < an_size; i++) 
    an_evol(gates[size+i], i);

}

void QSystem::evol(std::string u, size_t qbit) {
  auto size_m = log2(gate.cget(u).rows());

  for (size_t i = qbit; i < qbit+size_m; i++) {
    if ((i < size and ops[i] != 'I') or (i >= size and an_ops[size-i] != 'I')) {
      sync();
      break;
    } 
  }

  mops.push_back(u);
  char idx = mops.size()-1;
  
  for (size_t i = qbit; i < qbit+size_m; i++) {
    if (i < size) 
      ops[i] = idx;
    else 
      an_ops[size-i] = idx;
  }

  syncc = false;
}

void QSystem::sync() {
  sp_cx_mat evolm;
  
  size_t i;
  if (int(ops[0]) > 31) {
    evolm = gate.get(ops[0]);
    i = 1;
  } else {
    evolm = gate.cget(mops[int(ops[0])]);
    i = log2(evolm.rows());
  } 

  for (; i < size; i++) {
    if (int(ops[i]) > 31) {
      evolm = kron(evolm, gate.get(ops[i])); 
    } else {
      auto &m = gate.cget(mops[int(ops[i])]);
      evolm = kron(evolm, m);
      i += log2(m.rows())-1;
    }
  }

  for (i -= size; i < an_size; i++) {
    if (int(ops[i]) > 31) {
      evolm = kron(evolm, gate.get(an_ops[i])); 
    } else {
      auto &m = gate.cget(mops[i]);
      evolm = kron(evolm, m);
      i += log2(m.rows())-1;
    }
  }

  if (state == "pure")
    qbits = evolm*qbits;
  else if (state == "mix")
    qbits = evolm*qbits*evolm.adjoint();

  std::memset(ops, 'I', size*sizeof(char));
  std::memset(an_ops, 'I', an_size*sizeof(char));
  mops.clear();

  syncc = true;
}

void QSystem::clar() {
  std::memset(ops, 'I', size*sizeof(char));
  std::memset(an_ops, 'I', an_size*sizeof(char));
  mops.clear();

  syncc = true;
}

void QSystem::cnot(size_t target, vec_size control) {
  if (not syncc) sync();

  auto eyesize = 1l << (size+an_size);
  sp_cx_mat cnotm{eyesize, eyesize};

  for (size_t i = 0; i < (1lu << (size+an_size)); i++) {
    bool cond = true;

    for (size_t k = 0; (k < control.size()) and cond; k++)
      cond = cond and ((1ul << (size+an_size-control[k]-1)) & i);

    if (cond)
      cnotm.coeffRef(i, i ^ (1ul  << (size+an_size-target-1))) = 1; 
    else 
      cnotm.coeffRef(i, i) = 1;
  }

  if (state == "pure")
    qbits = cnotm*qbits;
  else if (state == "mix")
    qbits = cnotm*qbits*cnotm;
}

void QSystem::measure(size_t qbit) {
  if (qbit >= size+an_size)
    throw std::out_of_range("Argument 'qbit' must be in range  [0, "
        + std::to_string(size+an_size-1) + "].");

  if (not syncc) sync();

  double pm = 0;
  if (state == "pure") {
    for (it_mat i(qbits, 0); i; ++i) {
      if (~i.row() & 1ul << (size+an_size-qbit-1)) {
        auto valor = abs(i.value());
        pm += valor*valor;
      }
    }
  } else if (state == "mix") {
    sp_cx_vec m_aux = qbits.diagonal().sparseView(); 
    for (it_vec i(m_aux); i; ++i) 
      if (~i.row() & 1ul << (size+an_size-qbit-1)) 
        pm += i.value().real();
  }

  auto result = [&](Bit mea, double pm) {
    auto lnot = mea == zero? [](size_t i) { return ~i; } 
                           : [](size_t i) { return i; };
    if (qbit < size) bits[qbit] = mea;
      else an_bits[qbit-size] = mea;

    sp_cx_mat qbitsm{1u << (size+an_size), state == "pure"? 1 : 1u << (size+an_size)};

    if (state == "pure") {
      for (it_mat i(qbits, 0); i ; ++i) 
        if (lnot(i.row()) & 1ul << (size+an_size-qbit-1))  
          qbitsm.coeffRef(i.row(), 0) = i.value()/sqrt(pm);
    } else if (state == "mix") {
      for (auto k = 0l; k < qbits.outerSize(); ++k) {
        for (it_mat i(qbits, k); i; ++i) 
          if (lnot(i.row()) & 1ul << (size+an_size-qbit-1) 
              and lnot(i.col()) & 1ul << (size+an_size-qbit-1)) 
            qbitsm.coeffRef(i.row(), i.col()) = i.value()/pm;
      }
    }

    return qbitsm;
  };
  
  if (pm != 0 and double(std::rand()) / double(RAND_MAX) <= pm) 
    qbits = result(zero, pm);
  else 
    qbits = result(one, 1.0 - pm);
}

void QSystem::measure(size_t qbegin, size_t qend) {
  for (size_t i = qbegin; i < qend; i++) 
    measure(i);
}

void QSystem::measure_all() {
  for (size_t i = 0; i < size+an_size; i++) 
    measure(i);
}

