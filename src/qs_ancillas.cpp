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

void QSystem::add_ancillas(size_t an_num) {
  if (an_num == 0) 
    throw std::invalid_argument{"Argument \'an_num\' must be greater than 0"};
  if (not syncc) sync();

  an_size = an_num;
  an_ops = new char[an_size];
  an_bits = new Bit[an_size]();

  std::memset(an_ops, 'I', an_size*sizeof(char));

  sp_cx_mat an_qbits{1l << an_size, state == "mix" ? 1l << an_size : 1};
  an_qbits.coeffRef(0,0) = 1;
  qbits = kron(qbits, an_qbits);
}

void QSystem::rm_ancillas() {
  if (an_size == 0) 
    throw std::logic_error{"There are no ancillas on the system"};
  if (not syncc) sync();

  auto tr_pure = [&]() {
    auto sizet = 1l << (size+an_size-1);
    sp_cx_mat qbitst{sizet, 1};

    if (an_bits[an_size-1] == none) 
      an_measure(an_size-1);

    for (it_mat i(qbits, 0); i; ++i) 
      qbitst.coeffRef(i.row() >> 1, 0) += i.value();

    return qbitst;
  };

  auto tr_mix = [&]() {
    auto sizet = 1l << (size+an_size-1);
    sp_cx_mat qbitst{sizet, sizet};

    for (auto k = 0l; k < qbits.outerSize(); ++k) {
      for (it_mat i(qbits, k); i; ++i) 
        qbitst.coeffRef(i.row() >> 1, i.col() >> 1) += i.value();
    }

    return qbitst;
  };

  while (an_size) {
    if (state == "pure")
      qbits = tr_pure();
    else if (state == "mix")
      qbits = tr_mix();
    an_size--;
  }

  delete[] an_ops;
  an_ops = nullptr;
  delete[] an_bits;
  an_bits = nullptr;
}

void QSystem::an_evol(char gate, size_t qbit) {
  if (an_ops[qbit] != 'I') sync();
  
  an_ops[qbit] = gate;
  syncc = false;
}

void QSystem::an_evol(char gate, size_t qbegin, size_t qend) {
  for (size_t qbit = qbegin; qbit < qend; qbit++){
    if (an_ops[qbit] != 'I') sync();
    
    an_ops[qbit] = gate;
  }
  syncc = false;
}

void QSystem::an_measure(size_t qbit) {
  if (qbit >= an_size)
    throw std::out_of_range("Argument 'qbit' must be in range  [0, "
        + std::to_string(an_size-1) + "].");

  measure(size+qbit);
}

void QSystem::an_measure(size_t qbegin, size_t qend) {
  for (size_t qbit = qbegin; qbit < qend; qbit++) 
    measure(size+qbit);
}

