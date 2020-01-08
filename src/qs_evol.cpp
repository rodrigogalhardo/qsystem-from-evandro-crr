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
#include <algorithm>

using namespace arma;
using namespace std::complex_literals;

/******************************************************/
void QSystem::evol(char gate,
                 size_t qbit, 
                 size_t count,
                   bool invert) {
  valid_qbit("qbit", qbit);
  valid_count(qbit, count);
  if (_representation != "bitwise") {

    sync(qbit, qbit+count);
    for (size_t i = 0; i < count; i++) {
      ops(qbit+i).tag = Gate_aux::GATE_1;
      ops(qbit+i).data = gate;
      ops(qbit+i).invert = invert;
    }

    _sync = false;
  } else {
    for (size_t k = 0; k < count; k++) {
      switch (gate) {
        case 'I': return;
        case 'X': evol_x(qbit+k); break;
        case 'Y': evol_y(qbit+k); break;
        case 'Z': evol_z(qbit+k); break;
        case 'H': evol_h(qbit+k); break;
        case 'S': evol_s(qbit+k, invert); break;
        case 'T': evol_t(qbit+k, invert); break;
        }
    }
  }
}

/******************************************************/
void QSystem::rot(char axis, double angle, size_t qbit, size_t count) {
  valid_gate("axis", axis);
  valid_count(qbit, count);

  if (_representation != "bitwise") {
    sync(qbit, qbit+count);
    for (size_t i = 0; i < count; i++) {
      ops(qbit+i).tag = Gate_aux::R;
      ops(qbit+i).data = r_pair{axis, angle};
    }
    _sync = false;
  } else {
    for (size_t k = 0; k < count; k++) {
      dict bw_tmp;
      switch (axis) {
        case 'X':
          for (auto &i : bwqbits) {
            bw_tmp[i.first] += i.second*cos(angle/2);
            if (std::abs(bw_tmp[i.first]) < 1e-10) 
              bw_tmp.erase(i.first);
            
            size_t j = i.first ^ (1ul << (size()-qbit+k-1));
            bw_tmp[j] -= i.second*sin(angle/2)*complex(0, 1);
            if (std::abs(bw_tmp[j]) < 1e-10) 
              bw_tmp.erase(j);
          }
          bwqbits.swap(bw_tmp);

          break;

        case 'Y':
          for (auto &i : bwqbits) {
            bw_tmp[i.first] += i.second*cos(angle/2);
            if (std::abs(bw_tmp[i.first]) < 1e-10) 
              bw_tmp.erase(i.first);
            
            size_t j = i.first ^ (1ul << (size()-qbit+k-1));
            if (i.first & (1ul << (size()-qbit+k-1))) 
              bw_tmp[j] -= i.second*sin(angle/2);
            else 
              bw_tmp[j] += i.second*sin(angle/2);

            if (std::abs(bw_tmp[j]) < 1e-10) 
              bw_tmp.erase(j);
          }
          bwqbits.swap(bw_tmp);
 
          break;
        case 'Z':
          for (auto &i : bwqbits) 
            if (i.first & (1ul << (size()-qbit+k-1)))
              bwqbits[i.first] *= complex(cos(angle/2), sin(angle/2));
            else 
              bwqbits[i.first] *= -complex(cos(angle/2), sin(angle/2));
          break;
      }
    }
  }
}

/******************************************************/
void QSystem::u3(double theta, 
                 double phi,
                 double lambd,
                 size_t qbit, 
                 size_t count) {
  valid_count(qbit, count);
  if (_representation != "bitwise") {
    sync(qbit, qbit+count);
    for (size_t i = 0; i < count; i++) {
      ops(qbit+i).tag = Gate_aux::U3;
      ops(qbit+i).data = u3_tuple{theta, phi, lambd};
    }
    _sync = false;
  } else {
    for (size_t k = 0; k < count; k++) {
      dict bw_tmp;
      for (auto &i : bwqbits) {
        size_t j = i.first ^ (1ul << (size()-qbit+k-1));
        if (i.first & (1ul << (size()-qbit+k-1))) {
          bw_tmp[i.first] += i.second*complex(cos(lambd+phi)*cos(theta/2), 
                                              cos(theta/2)*sin(lambd+phi));
          if (std::abs(bw_tmp[i.first]) < 1e-10) 
            bw_tmp.erase(i.first);

          bw_tmp[j] += i.second*(-complex(cos(lambd)*sin(theta/2), 
                                          sin(lambd)*sin(theta/2)));
          if (std::abs(bw_tmp[j]) < 1e-10) 
            bw_tmp.erase(j);
        } else {
          bw_tmp[i.first] += i.second*cos(theta/2); 
          if (std::abs(bw_tmp[i.first]) < 1e-10) 
            bw_tmp.erase(i.first);

          bw_tmp[j] += i.second*complex(cos(phi)*sin(theta/2), 
                                        sin(phi)*sin(theta/2));
          if (std::abs(bw_tmp[j]) < 1e-10) 
            bw_tmp.erase(j);
        }
      }
      bwqbits.swap(bw_tmp);
    }
  }
}

/******************************************************/
void QSystem::u2(double phi,
                 double lambd,
                 size_t qbit, 
                 size_t count) {
  valid_count(qbit, count);
  if (_representation != "bitwise") {
    sync(qbit, qbit+count);
    for (size_t i = 0; i < count; i++) {
      ops(qbit+i).tag = Gate_aux::U3;
      ops(qbit+i).data = u3_tuple{acos(-1)/2, phi, lambd};
    }
    _sync = false;
  } else {
    for (size_t k = 0; k < count; k++) {
      dict bw_tmp;
      for (auto &i : bwqbits) {
        size_t j = i.first ^ (1ul << (size()-qbit+k-1));
        if (i.first & (1ul << (size()-qbit+k-1))) {
          bw_tmp[i.first] += i.second*(complex(cos(lambd+phi), 
                                              sin(lambd+phi))/sqrt(2));
          if (std::abs(bw_tmp[i.first]) < 1e-10) 
            bw_tmp.erase(i.first);

          bw_tmp[j] += i.second*(-complex(cos(lambd), 
                                          sin(lambd))/sqrt(2));
          if (std::abs(bw_tmp[j]) < 1e-10) 
            bw_tmp.erase(j);
        } else {
          bw_tmp[i.first] += i.second/sqrt(2); 
          if (std::abs(bw_tmp[i.first]) < 1e-10) 
            bw_tmp.erase(i.first);

          bw_tmp[j] += i.second*(complex(cos(phi), 
                                        sin(phi))/sqrt(2));
          if (std::abs(bw_tmp[j]) < 1e-10) 
            bw_tmp.erase(j);
        }
      }
      bwqbits.swap(bw_tmp);
    }
  }
}

/******************************************************/
void QSystem::u1(double lambd,
                 size_t qbit, 
                 size_t count) {
  valid_count(qbit, count);
  if (_representation != "bitwise") {
    sync(qbit, qbit+count);
    for (size_t i = 0; i < count; i++) {
      ops(qbit+i).tag = Gate_aux::U3;
      ops(qbit+i).data = u3_tuple{0, 0, lambd};
    }
    _sync = false;
  } else {
    for (size_t k = 0; k < count; k++) {
      for (auto &i : bwqbits) 
        if (i.first & (1ul << (size()-qbit+k-1)))
          bwqbits[i.first] *= complex(cos(lambd), sin(lambd));
    }
  }
}

/******************************************************/
void QSystem::apply(Gate gate, size_t qbit, size_t count, bool invert) {
  size_t size_n = log2(gate.get_mat()->n_rows);
  valid_count(qbit, count, size_n);
  if (_representation != "bitwise") {
    sync(qbit, qbit+count*size_n);
    for (size_t i = 0; i < count; i++) {
      size_t index = qbit+i*size_n;
      fill(Gate_aux::GATE_N, index, size_n);
      ops(index).data = gate.get_mat();
      ops(index).invert = invert;
    }
  } else {
    dict bw_tmp;
    for (auto &i : bwqbits) {
      // i.first = x|y|z
      size_t x = i.first & (((1ul << qbit)-1) << (size()-qbit));
      size_t y = i.first >> (size()-qbit-size_n);
      y = y & ((1ul << size_n)-1);
      size_t z = i.first & ((1ul << (size()-qbit-size_n))-1);
      auto &setu = gate.get_bwgate(y);
      for (auto &j : setu) {
        size_t xjz = x|(j.second << (size()-qbit-size_n))|z;
        bw_tmp[xjz] += i.second*j.first;
        if (std::abs(bw_tmp[xjz]) < 1e-10) 
          bw_tmp.erase(xjz);
      }
    }
    bwqbits.swap(bw_tmp);
  }
}

/******************************************************/
void QSystem::cnot(size_t target, vec_size_t control) {
  valid_qbit("target", target);
  valid_control(control);
  if (_representation != "bitwise") {
    auto [size_n, minq] = cut(target, control);
    fill(Gate_aux::CNOT, minq, size_n);
    ops(minq).data = cnot_pair{target, control};
  } else {
    dict bw_tmp;
    for (auto &i : bwqbits) {
      bool tmp_ctrl = true;
      for (auto ctrl : control) {
        tmp_ctrl = tmp_ctrl and (i.first & (1ul << (size()-ctrl-1)));
      }
      if (tmp_ctrl) {
        size_t j = i.first ^ (1ul << (size()-target-1));
        bw_tmp[j] = i.second;
      } else {
        bw_tmp[i.first] = i.second;
      }
    }
    bwqbits.swap(bw_tmp);
  }
}

/******************************************************/
void QSystem::cphase(complex phase, size_t target, vec_size_t control) {
  valid_qbit("target", target);
  valid_phase(phase);
  valid_control(control);

  if (_representation != "bitwise") {
    auto [size_n, minq] = cut(target, control);
    fill(Gate_aux::CPHASE, minq, size_n);
    ops(minq).data = cph_tuple{phase, target, control};
  } else {
    dict bw_tmp;
    for (auto &i : bwqbits) {
      bool tmp_ctrl = true;
      for (auto ctrl : control) {
        tmp_ctrl = tmp_ctrl and (i.first & (1ul << (size()-ctrl-1)));
      }
      if (tmp_ctrl) {
        size_t j = i.first ^ (1ul << (size()-target-1));
        bw_tmp[j] = i.second*phase;
      } else {
        bw_tmp[i.first] = i.second;
      }
    }
    bwqbits.swap(bw_tmp);
  }
}

/******************************************************/
void QSystem::swap(size_t qbit_a, size_t qbit_b) {
  valid_swap(qbit_a, qbit_b);
  
  if (_representation != "bitwise") {
    if (qbit_a == qbit_b) return;
    size_t a = qbit_a < qbit_b? qbit_a :  qbit_b;
    size_t b = qbit_a > qbit_b? qbit_a :  qbit_b;
    fill(Gate_aux::SWAP, a, b-a+1);
  } else {
    dict bw_tmp;
    for (auto &i : bwqbits) {
      bool bit_a = i.first & (1ul << (size()-qbit_a-1));
      bool bit_b = i.first & (1ul << (size()-qbit_b-1));
      if (bit_a != bit_b) {
        size_t j = i.first ^ (1ul << (size()-qbit_a-1));
        j = j ^ (1ul << (size()-qbit_b-1));
        bw_tmp[j] = i.second;
      } else {
        bw_tmp[i.first] = i.second;
      }
    }
    bwqbits.swap(bw_tmp);
  }
}

/******************************************************/
void QSystem::qft(size_t qbit_begin, size_t qbit_end, bool invert) {
  valid_range(qbit_begin, qbit_end);
  auto size_n = qbit_end-qbit_begin;
  if (representation() != "bitwise") {
    fill(Gate_aux::QFT, qbit_begin, qbit_end-qbit_begin);
    ops(qbit_begin).invert = invert;
  } else {
    double pi = acos(-1);
    complex w = std::exp((2*pi*1i)/(pow(2, size_n)));
    dict bw_tmp;
   
    for (auto &i : bwqbits) {
      // i.first = x|y|z
      size_t x = i.first & (((1ul << qbit_begin)-1) << (size()-qbit_begin));
      size_t y = i.first >> (size()-qbit_begin-size_n);
      y = y & ((1ul << size_n)-1);
      size_t z = i.first & ((1ul << (size()-qbit_begin-size_n))-1);
      
      for (size_t j = 0; j < (1ul << size_n); j++) {
        auto val = (1/sqrt(1 << size_n))*pow(w, y*j);
        auto xjz = x|(j << (size()-qbit_begin-size_n))|z;
        bw_tmp[xjz] += val*i.second;
        if (std::abs(bw_tmp[xjz]) < 1e-10) 
          bw_tmp.erase(xjz);
      }
    }

    bwqbits.swap(bw_tmp);
  }
}

/******************************************************/
void QSystem::sync() {
  if (_sync) return;

  sp_cx_mat evolm;
  
  evolm = get_gate(ops(0));

  for (size_t i = ops(0).size; i < size(); i += ops(i).size) {
    evolm = kron(evolm, get_gate(ops(i)));
  }

  if (_representation == "vector")
    qbits = evolm*qbits;
  else if (_representation == "matrix")
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
void QSystem::sync(size_t qbit_begin, size_t qbit_end) {
  for (size_t i = qbit_begin; i < qbit_end; i++) {
    if (ops(i).busy()) {
      sync();
      break;
    } 
  }
}

/******************************************************/
QSystem::Gate_aux& QSystem::ops(size_t index) {
  return index < _size? _ops[index] : an_ops[index-_size];
}


/******************************************************/
sp_cx_mat QSystem::get_gate(Gate_aux &op) {
  auto get = [&]() {
      switch (op.tag) {
      case Gate_aux::GATE_1:
        return gates[std::get<char>(op.data)];
      case Gate_aux::GATE_N:
        return *std::get<mat_ptr>(op.data);
      case Gate_aux::R:
        return make_rot(std::get<r_pair>(op.data).first,
                      std::get<r_pair>(op.data).second);
      case Gate_aux::U3:
        return make_u3(std::get<0>(std::get<u3_tuple>(op.data)),
                       std::get<1>(std::get<u3_tuple>(op.data)),
                       std::get<2>(std::get<u3_tuple>(op.data)));
      case Gate_aux::CNOT:
        return make_cnot(std::get<cnot_pair>(op.data).first,
                         std::get<cnot_pair>(op.data).second,
                         op.size);
      case Gate_aux::CPHASE:
        return make_cphase(std::get<0>(std::get<cph_tuple>(op.data)),
                           std::get<1>(std::get<cph_tuple>(op.data)),
                           std::get<2>(std::get<cph_tuple>(op.data)),
                           op.size);
      case Gate_aux::SWAP:
        return make_swap(op.size);
      default:
        return make_qft(op.size);
      }
  };

  if (op.invert) 
    return get().t();
  else 
    return get();
}

/******************************************************/
cut_pair QSystem::cut(size_t &target, vec_size_t &control) {
  size_t maxq = std::max(target,
                         *std::max_element(control.begin(),
                                           control.end()));
  size_t minq = std::min(target,
                         *std::min_element(control.begin(),
                                           control.end()));
  size_t size_n = maxq-minq+1;
  for (auto &i : control) 
    i -= minq;
  target -= minq;
  return std::make_pair(size_n, minq);
}

/******************************************************/
void QSystem::fill(Gate_aux::Tag tag, size_t qbit, size_t size_n) {
  sync(qbit, qbit+size_n);

  ops(qbit).tag = tag;
  ops(qbit).size = size_n;
 
  for (size_t i = qbit+1; i < qbit+size_n; i++)
    ops(i).tag = tag;

  _sync = false;
}

/******************************************************/
void QSystem::evol_h(size_t qbit) {
  dict bw_tmp;

  for (auto &i : bwqbits) {
    if (i.first & (1ul << (size()-qbit-1)))
      bw_tmp[i.first] -= i.second/std::sqrt(2);
    else 
      bw_tmp[i.first] += i.second/std::sqrt(2);
    if (std::abs(bw_tmp[i.first]) < 1e-10) 
      bw_tmp.erase(i.first);

    size_t j = i.first ^ (1ul << (size()-qbit-1));
    bw_tmp[j] += i.second/std::sqrt(2);
    if (std::abs(bw_tmp[j]) < 1e-10) 
      bw_tmp.erase(j);
  }

  bwqbits.swap(bw_tmp);
}

/******************************************************/
void QSystem::evol_x(size_t qbit) {
  dict bw_tmp;
  for (auto &i : bwqbits) {
    size_t j = i.first ^ (1ul << ((size()-qbit-1)));
    bw_tmp[j] = i.second; 
  }

  bwqbits.swap(bw_tmp);
}

/******************************************************/
void QSystem::evol_y(size_t qbit) {
  dict bw_tmp;
  for (auto& i : bwqbits) {
    size_t j = i.first ^ (1ul << (size()-qbit-1));
    if (i.first & (1ul << (size()-qbit-1)))
      bw_tmp[j] = i.second*complex(0, -1);
    else
      bw_tmp[j] = i.second*complex(0, 1);
  }

  bwqbits.swap(bw_tmp);
}

/******************************************************/
void QSystem::evol_z(size_t qbit) {
  for (auto &i : bwqbits) 
    if (i.first & (1ul << (size()-qbit-1)))
      bwqbits[i.first] *= -1;
}

/******************************************************/
void QSystem::evol_s(size_t qbit, bool invert) {
  for (auto &i : bwqbits) {
    if (i.first & (1ul << (size()-qbit-1))) {
      if (invert)
        bwqbits[i.first] *= complex(0, -1);
      else
        bwqbits[i.first] *= complex(0, 1);
    }
  }
}

/******************************************************/
void QSystem::evol_t(size_t qbit, bool invert) {
  for (auto &i : bwqbits) {
    if (i.first & (1ul << (size()-qbit-1))) {
      if (invert)
        bwqbits[i.first] *= complex(1/std::sqrt(2), -1/std::sqrt(2));
      else 
        bwqbits[i.first] *= complex(1/std::sqrt(2), 1/std::sqrt(2));
    }
  }
}
