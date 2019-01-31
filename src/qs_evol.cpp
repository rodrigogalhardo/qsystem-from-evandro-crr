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

using namespace arma;

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
  auto size_m = log2(gate.cget(u).n_rows);

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
    i = log2(evolm.n_rows);
  } 

  for (; i < size; i++) {
    if (int(ops[i]) > 31) {
      evolm = kron(evolm, gate.get(ops[i])); 
    } else {
      auto &m = gate.cget(mops[int(ops[i])]);
      evolm = kron(evolm, m);
      i += log2(m.n_rows)-1;
    }
  }

  for (i -= size; i < an_size; i++) {
    if (int(ops[i]) > 31) {
      evolm = kron(evolm, gate.get(an_ops[i])); 
    } else {
      auto &m = gate.cget(mops[i]);
      evolm = kron(evolm, m);
      i += log2(m.n_rows)-1;
    }
  }

  if (state == "pure")
    qbits = evolm*qbits;
  else if (state == "mix")
    qbits = evolm*qbits*evolm.t();

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

  size_t eyesize = 1ul << (size+an_size);
  sp_cx_mat cnotm{eyesize, eyesize};

  for (size_t i = 0; i < (1lu << (size+an_size)); i++) {
    bool cond = true;

    for (size_t k = 0; (k < control.size()) and cond; k++)
      cond = cond and ((1ul << (size+an_size-control[k]-1)) & i);

    if (cond)
      cnotm(i, i ^ (1ul  << (size+an_size-target-1))) = 1; 
    else 
      cnotm(i, i) = 1;
  }

  if (state == "pure")
    qbits = cnotm*qbits;
  else if (state == "mix")
    qbits = cnotm*qbits*cnotm;
}

/*
void QSystem::ctrl(std::string gates, vec_size control) {
  auto par = [&](size_t x) {
    for (int i = 32; i > 0; i /= 2)
      x ^= x >> i;
    return x & 1; 
  };
  if (not syncc) sync();

  size_t eyesize = 1ul << (size+an_size);
  sp_cx_mat nqbits{eyesize, eyesize};

  size_t x = 0;
  size_t z = 0;
  for (size_t i = 0; i < size+an_size; i++) {
    if (gates[i] == 'X')
      x |= 1ul << ((size+an_size)-i-1);
    if (gates[i] == 'Z')
      z |= 1ul << ((size+an_size)-i-1);
  }

  auto cond = [&](size_t ij) {
    bool c = true;
    for (size_t k = 0; (k < control.size()) and c; k++)
      c = c and ((1ul << (size+an_size-control[k]-1)) & ij);
    if (c) 
      return std::make_pair(ij^x, pow(-1, par(ij & z)));
    else 
      return std::make_pair(ij, 1.0);
  };

  for (auto i = qbits.begin(); i != qbits.end(); ++i) {
    auto [ii, fi] = cond(i.row());
    auto [ij, fj] = cond(i.col());
    
    nqbits(ii, ij) = ((std::complex<double>)*i)*fi*fj;
  }

  qbits = nqbits;
}
*/

void QSystem::measure(size_t qbit) {
  if (qbit >= size+an_size)
    throw std::out_of_range("Argument 'qbit' must be in range  [0, "
        + std::to_string(size+an_size-1) + "].");

  if (not syncc) sync();

  double pm = 0;
  if (state == "pure") {
    for (auto i = qbits.begin(); i != qbits.end(); ++i) {
      if (~i.row() & 1ul << (size+an_size-qbit-1)) {
        auto valor = abs((cx_double) *i);
        pm += valor*valor;
      }
    }
  } else if (state == "mix") {
    sp_cx_mat m_aux = qbits.diag(); 
    for (auto i = m_aux.begin(); i != m_aux.end(); ++i) 
      if (~i.row() & 1ul << (size+an_size-qbit-1)) 
        pm += (*i).real();
  }

  auto result = [&](Bit mea, double pm) {
    auto lnot = mea == zero? [](size_t i) { return ~i; } 
                           : [](size_t i) { return i; };
    if (qbit < size) bits[qbit] = mea;
      else an_bits[qbit-size] = mea;

    sp_cx_mat qbitsm{1ul << (size+an_size), state == "pure"? 1 : 1ul << (size+an_size)};

    if (state == "pure") {
      for (auto i = qbits.begin(); i != qbits.end(); ++i) 
        if (lnot(i.row()) & 1ul << (size+an_size-qbit-1))  
          qbitsm(i.row(), 0) = (cx_double)(*i)/sqrt(pm);
    } else if (state == "mix") {
      for (auto i = qbits.begin(); i != qbits.end(); ++i) 
        if (lnot(i.row()) & 1ul << (size+an_size-qbit-1) 
            and lnot(i.col()) & 1ul << (size+an_size-qbit-1)) 
          qbitsm(i.row(), i.col()) = (cx_double)(*i)/pm;
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

