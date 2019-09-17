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
void QSystem::add_ancillas(size_t num_qbits, size_t init) {
  valid_init(init, num_qbits);
  if (num_qbits == 0) {
    throw std::invalid_argument{"\'an_num\' argument must be greater than 0"};
  } else if (an_size != 0) {
    sstr err;
    err << "There are already ancillas in the system, you can not add more";
    throw std::invalid_argument{err.str()};
  }

  sync();

  an_size = num_qbits;

  if (representation() != "bitwise") {
    an_ops = new Gate_aux[an_size]();

    sp_cx_mat an_qbits{1ul << an_size, _representation == "matrix" ? 1lu << an_size : 1};
    an_qbits(init, _representation == "matrix" ? init : 0) = 1;
    qbits = kron(qbits, an_qbits);
  } else {
    dict bw_tmp;
    for (auto &i : bwqbits) {
      auto j = (i.first << num_qbits) | init;
      bw_tmp[j] = i.second;
    }
    bwqbits.swap(bw_tmp);
  }

  an_bits = new Bit[an_size]();
}

/******************************************************/
void QSystem::rm_ancillas() {
  if (an_size == 0) 
    throw std::logic_error{"There are no ancillas on the system"};
  sync();

  auto tr_pure = [&]() {
    auto sizet = 1ul << (size()-1);
    sp_cx_mat qbitst{sizet, 1};

    if (an_bits[an_size-1] == NONE) 
      measure(_size+an_size-1);

    for (auto i = qbits.begin(); i != qbits.end(); ++i) 
      qbitst(i.row() >> 1, 0) += *i;

    return qbitst;
  };

  auto tr_mix = [&]() {
    auto sizet = 1ul << (size()-1);
    sp_cx_mat qbitst{sizet, sizet};

    for (auto i = qbits.begin(); i != qbits.end(); ++i) {
      if ((i.row() % 2) == (i.col() % 2))
        qbitst(i.row() >> 1, i.col() >> 1) += *i;
    }

    return qbitst;
  };

  if (representation() != "bitwise") {
    while (an_size) {
      if (_representation == "vector")
        qbits = tr_pure();
      else if (_representation == "matrix")
        qbits = tr_mix();
      an_size--;
    }

    delete[] an_ops;
    an_ops = nullptr;
  } else {
    measure(_size, an_size);
    dict bw_tmp;
    
    for (auto &i : bwqbits) {
      auto j = i.first >> an_size;
      bw_tmp[j] = i.second;
    }

    bwqbits.swap(bw_tmp);
    an_size = 0;
  }

  delete[] an_bits;
  an_bits = nullptr;
}

