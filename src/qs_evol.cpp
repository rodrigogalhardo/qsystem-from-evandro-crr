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
void QSystem::evol(char gate,
                 size_t qbit, 
                 size_t count,
                   bool inver) {
  valid_qbit("qbit", qbit);

  valid_count(qbit, count);
  sync(qbit, qbit+count);
  for (size_t i = 0; i < count; i++) {
    ops(qbit+i).tag = Gate_aux::GATE_1;
    ops(qbit+i).data = gate;
    ops(qbit+i).inver = inver;
  }

  _sync = false;
}

/******************************************************/
void QSystem::r(char axis, double angle, size_t qbit, size_t count) {
  valid_gate("axis", axis);

  sync(qbit, qbit+count);
  for (size_t i = 0; i < count; i++) {
    ops(qbit+i).tag = Gate_aux::R;
    ops(qbit+i).data = r_pair{axis, angle};
  }
  _sync = false;
}

/******************************************************/
void QSystem::u3(double theta, 
                 double phi,
                 double lambd,
                 size_t qbit, 
                 size_t count) {
  sync(qbit, qbit+count);
  for (size_t i = 0; i < count; i++) {
    ops(qbit+i).tag = Gate_aux::U3;
    ops(qbit+i).data = u3_tuple{theta, phi, lambd};
  }
  _sync = false;
}

/******************************************************/
void QSystem::u2(double phi,
                 double lambd,
                 size_t qbit, 
                 size_t count) {
  sync(qbit, qbit+count);
  for (size_t i = 0; i < count; i++) {
    ops(qbit+i).tag = Gate_aux::U3;
    ops(qbit+i).data = u3_tuple{M_PI/2, phi, lambd};
  }
  _sync = false;
}

/******************************************************/
void QSystem::u1(double lambd,
                 size_t qbit, 
                 size_t count) {
  sync(qbit, qbit+count);
  for (size_t i = 0; i < count; i++) {
    ops(qbit+i).tag = Gate_aux::U3;
    ops(qbit+i).data = u3_tuple{0, 0, lambd};
  }
  _sync = false;
}

/******************************************************/
void QSystem::apply(Gate gate, size_t qbit, size_t count, bool inver) {
  auto size_n = log2(gate.get_mat()->n_rows);
  valid_count(qbit, count, size_n);
  sync(qbit, qbit+count*size_n);
  for (size_t i = 0; i < count; i++) {
    size_t index = qbit+i*size_n;
    fill(Gate_aux::GATE_N, index, size_n);
    ops(index).data = gate.get_mat();
    ops(index).inver = inver;
  }
}

/******************************************************/
void QSystem::cnot(size_t target, vec_size_t control) {
  valid_qbit("target", target);
  valid_control(control);

  auto [size_n, minq] = cut(target, control);
  fill(Gate_aux::CNOT, minq, size_n);
  ops(minq).data = cnot_pair{target, control};
}

/******************************************************/
void QSystem::cphase(complex phase, size_t target, vec_size_t control) {
  valid_qbit("target", target);
  valid_phase(phase);
  valid_control(control);

  auto [size_n, minq] = cut(target, control);
  fill(Gate_aux::CPHASE, minq, size_n);
  ops(minq).data = cph_tuple{phase, target, control};
}

/******************************************************/
void QSystem::swap(size_t qbit_a, size_t qbit_b) {
  valid_swap(qbit_a, qbit_b);

  if (qbit_a == qbit_b) return;
  size_t a = qbit_a < qbit_b? qbit_a :  qbit_b;
  size_t b = qbit_a > qbit_b? qbit_a :  qbit_b;
  fill(Gate_aux::SWAP, a, b-a+1);
}

/******************************************************/
void QSystem::qft(size_t qbegin, size_t qend, bool inver) {
  valid_range(qbegin, qend);

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
        return make_r(std::get<r_pair>(op.data).first,
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

  if (op.inver) 
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

