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

    sp_cx_mat qbitsm{1ul << (size+an_size),
                     state == "pure"? 1 : 1ul << (size+an_size)};

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

/******************************************************/
void QSystem::measure(size_t qbegin, size_t qend) {
  for (size_t i = qbegin; i < qend; i++) 
    measure(i);
}

/******************************************************/
void QSystem::measure_all() {
  for (size_t i = 0; i < size+an_size; i++) 
    measure(i);
}

