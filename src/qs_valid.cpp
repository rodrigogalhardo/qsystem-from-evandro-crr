/* MIT License
 * 
 * Copyright (c) 2019 Bruno Gouvêa Taketani <b.taketani@ufsc.br>
 * Copyright (c) 2019 Evandro Chagas Ribeiro da Rosa <ev.crr97@gmail.com>
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

/******************************************************/
void QSystem::valid_qbit(std::string name, size_t qbit) {
  if (qbit >= size()) {
      sstr err;
      err << "\'" << name << "\' argument should be in the range of 0 to "
          << (size()-1);
      throw std::invalid_argument{err.str()};
  }
}

/******************************************************/
void QSystem::valid_count(size_t qbit, size_t count, size_t size_n) {
  if (count == 0 and qbit+count*size_n <= size()) {
      sstr err;
      err << "\'cout\' argument should be greater than 0 "
          << "and \'qbit+count\' suld be in the range of 0 to "
          << size();
      throw std::invalid_argument{err.str()};
  }
}

/******************************************************/
void QSystem::valid_control(vec_size_t &control) {
  if (control.size() == 0) {
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
}

/******************************************************/
void QSystem::valid_phase(complex phase) {
  if (std::abs(std::abs(phase) - 1.0) > 1e-14) {
    sstr err;
    err << "abs(phase) must be equal to 1";
    throw std::invalid_argument{err.str()};
  }
}

/******************************************************/
void QSystem::valid_swap(size_t qbit_a, size_t qbit_b) {
  if (qbit_a >= size() or qbit_b >= size()) {
      sstr err;
      err << "Arguments \'qbit_a\' and \'qbit_b\' should be in the "
          << "range of 0 to " << (size()-1);
      throw std::invalid_argument{err.str()};
  }
}

/******************************************************/
void QSystem::valid_range(size_t qbegin, size_t qend) {
  if (qbegin >= size() or qend > size() or qbegin <= qend) {
      sstr err;
      err << "\'qbegin\' argument should be in the "
          << "range of 0 to " << (size()-1)
          << " and argument \'qend\' should be greater than 0 "
          << "and in the range of 0 to " << size();
      throw std::invalid_argument{err.str()};
  }
}

/******************************************************/
void QSystem::valid_gate(char gate) {
  if (not (gate == 'X' or gate == 'Y' or gate == 'Z')) {
    sstr err;
    err << "\'gate\' argument must be equal to \'X\', \'Y\' or \'X\'";
    throw std::invalid_argument{err.str()};    
  }
}

/******************************************************/
void QSystem::valid_p(double p) {
  if (p < 0 or p > 1) {
    sstr err;
    err << "\'p\' argument should be in the range of 0.0 to 1.0";
    throw std::invalid_argument{err.str()};
  }
}

void QSystem::valid_state() {
  if (_state == "vector") {
    sstr err;
    err << "\'state\' must be in \"matrix\" to apply this channel";
    throw std::runtime_error{err.str()};
  }
}

void QSystem::valid_krau(vec_str &kraus) {
  size_t ksize = kraus[0].size();
  for (auto& k : kraus) {
    if (k.size() != ksize) {
      sstr err;
      err << "All \'kraus\' operators must have the same size";
      throw std::runtime_error{err.str()};
    }
  }
}
