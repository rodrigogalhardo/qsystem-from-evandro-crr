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
void QSystem::evol(char gate, size_t qbit) {
  if (ops[qbit].tag != Op::NONE) sync();

  ops[qbit].tag = Op::GATE_1;
  ops[qbit].data = gate;
  syncc = false;
}

/******************************************************/
void QSystem::evol(char gate, size_t qbegin, size_t qend) {
  for (size_t i = qbegin; i < qend; i++) 
    evol(gate, i);
}

/******************************************************/
void QSystem::evol(std::string gates) {
  if (gates.size() != size+an_size) 
    throw std::length_error{"String 'gates' must be " +
        std::to_string(size+an_size) + " characters long"};
  
  for (size_t i = 0; i < size; i++) 
    evol(gates[i], i);

  for (size_t i = 0; i < an_size; i++) 
    an_evol(gates[size+i], i);
}

/******************************************************/
void QSystem::evol(std::string u, size_t qbit) {
  auto size_n = log2(gate.cget(u).n_rows);

  sync(qbit, qbit+size_n);

  ops[qbit].tag = Op::GATE_N;
  ops[qbit].data = u;
  ops[qbit].size = size_n;
  
  syncc = false;
}

/******************************************************/
void QSystem::cnot(size_t target, vec_size control) {
  auto [size_n, minq] = cut(target, control);

  sync(minq, minq+size_n);

  ops[minq].tag = Op::CNOT;
  ops[minq].data = cnot_pair{target, control};
  ops[minq].size = size_n;

  syncc = false;
}

/******************************************************/
void QSystem::cphase(cx phase, size_t target, vec_size control) {
  auto [size_n, minq] = cut(target, control);

  sync(minq, minq+size_n);

  ops[minq].tag = Op::CPHASE;
  ops[minq].data = cph_tuple{phase, target, control};
  ops[minq].size = size_n;

  syncc = false;
}

/******************************************************/
void QSystem::swap(size_t qbit_a, size_t qbit_b) {
  if (qbit_a == qbit_b) return;

  size_t a = qbit_a < qbit_b? qbit_a :  qbit_b;
  size_t b = qbit_a > qbit_b? qbit_a :  qbit_b;

  size_t size_n = b-a+1;

  sync(a, a+size_n);
 
  ops[a].tag = Op::SWAP;
  ops[a].size = size_n;

  syncc = false;
}

/******************************************************/
void QSystem::qft(size_t qbegin, size_t qend) {
  sync(qbegin, qend);

  ops[qbegin].tag = Op::QFT;
  ops[qbegin].size = qend-qbegin;

  syncc = false;
}

/******************************************************/
void QSystem::sync() {
  sp_cx_mat evolm;
  
  evolm = get_gate(ops[0]);
  size_t i = ops[0].size;
  
  for (; i < size; i += ops[i].size) 
    evolm = kron(evolm, get_gate(ops[i]));

  for (i -= size; i < an_size; i += an_ops[i].size) 
    evolm = kron(evolm, get_gate(an_ops[i]));

  if (state == "pure")
    qbits = evolm*qbits;
  else if (state == "mix")
    qbits = evolm*qbits*evolm.t();

  ops = new (ops) Op[size];
  if (an_ops)
    an_ops = new (an_ops) Op[an_size];

  syncc = true;
}

/******************************************************/
void QSystem::sync(size_t qbegin, size_t qend) {
  for (size_t i = qbegin; i < qend; i++) {
    if ((i < size and ops[i].tag != Op::NONE) or
        (i >= size and an_ops[size-i].tag != Op::NONE)) {
      sync();
      break;
    } 
  }
}

/******************************************************/
void QSystem::clar() {
  ops = new (ops) Op[size];
  if (an_ops)
    an_ops = new (an_ops) Op[an_size];

  syncc = true;
}

