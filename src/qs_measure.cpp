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
void QSystem::measure(size_t qbit, size_t count) {
  valid_qbit("qibt", qbit);
  valid_count(qbit, count);

  sync();
  count += qbit;
  for (; qbit < count; qbit++) {
    double pm = 0;
    if (_representation == "vector") {
      for (auto i = qbits.begin(); i != qbits.end(); ++i) {
        if (~i.row() & 1ul << (size()-qbit-1)) {
          auto valor = abs((cx_double) *i);
          pm += valor*valor;
        }
      }
    } else if (_representation == "matrix") {
      sp_cx_mat m_aux = qbits.diag(); 
      for (auto i = m_aux.begin(); i != m_aux.end(); ++i) 
        if (~i.row() & 1ul << (size()-qbit-1)) 
          pm += (*i).real();
    } else {
      for (auto &i : bwqbits) {
        if (~i.first & 1ul << (size()-qbit-1)) {
          auto valor = abs(i.second);
          pm += valor*valor;
        }
      }
    }

    auto mea = pm != 0 and double(std::rand()) / double(RAND_MAX) <= pm? ZERO : ONE;
    pm = mea == ZERO? pm : 1.0 - pm;

    auto lnot = mea == ZERO? [](size_t i) { return ~i; } 
                            : [](size_t i) { return i; };

    if (qbit < _size) _bits[qbit] = mea;
      else an_bits[qbit-_size] = mea;

    if (representation() != "bitwise") {
      sp_cx_mat qbitsm{1ul << size(),
                       _representation == "vector"? 1 : 1ul << size()};

      if (_representation == "vector") {
        for (auto i = qbits.begin(); i != qbits.end(); ++i) 
          if (lnot(i.row()) & 1ul << (size()-qbit-1))  
            qbitsm(i.row(), 0) = (complex)(*i)/sqrt(pm);
      } else if (_representation == "matrix") {
        for (auto i = qbits.begin(); i != qbits.end(); ++i) 
          if (lnot(i.row()) & 1ul << (size()-qbit-1)
              and lnot(i.col()) & 1ul << (size()-qbit-1)) 
            qbitsm(i.row(), i.col()) = (complex)(*i)/pm;
      }

      qbits = qbitsm;
    } else {
      dict bw_tmp;
      for (auto &i : bwqbits) {
        if (lnot(i.first) & 1ul << (size()-qbit-1))  
          bw_tmp[i.first] = i.second/sqrt(pm);
      }
      bwqbits.swap(bw_tmp);
    }
  }
}

/******************************************************/
void QSystem::measure_all() {
  measure(0, size());
}

/******************************************************/
std::vector<int> QSystem::bits() {
  std::vector<int> vec;
  for (size_t i = 0; i < _size; i++)
    vec.push_back(_bits[i]);
  for (size_t i = 0; i < an_size; i++)
    vec.push_back(an_bits[i]);
  return vec;
}
