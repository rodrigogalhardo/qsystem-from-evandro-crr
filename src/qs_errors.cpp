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

void QSystem::flip(char gate, size_t qbit, double p) {
  if (not syncc) sync();

  if (state == "pure") {
    if (auto prand = double(std::rand()) / double(RAND_MAX); p != 0 and prand <= p) 
      evol(gate, qbit);
  } else if (state == "mix") {
    sp_cx_mat e0 = make_gate(this->gate.get(gate), qbit)*sqrt(p);
    size_t eyesize = 1ul << (size+an_size);
    sp_cx_mat e1 = eye<sp_cx_mat>(eyesize, eyesize)*sqrt(1.f-p);
    qbits = e0*qbits*e0 + e1*qbits*e1;
  }

}

void QSystem::amp_damping(size_t qbit, double p) {
  if (state == "pure")
    throw std::runtime_error{"\'state\' must be in \"mix\" to apply this channel."};

  if (not syncc) sync();

  sp_cx_mat e0 = make_gate(sp_cx_mat{cx_mat{{{{1, 0}, {0, 0}},
                                             {{0, 0}, {sqrt(1-p), 0}}}}}, qbit);
  sp_cx_mat e1 = make_gate(sp_cx_mat{cx_mat{{{{0, 0}, {sqrt(p), 0}},
                                             {{0, 0}, {0, 0}}}}}, qbit);
  qbits = e0*qbits*e0 + e1*qbits*e1;
}

void QSystem::dpl_channel(size_t qbit, double p) {
  if (state == "pure")
    throw std::runtime_error{"\'state\' must be in \"mix\" to apply this channel."};

  if (not syncc) sync();

  sp_cx_mat x = make_gate(gate.get('X'), qbit);
  sp_cx_mat y = make_gate(gate.get('Y'), qbit);
  sp_cx_mat z = make_gate(gate.get('Z'), qbit);
  
  qbits = (1-p)*qbits+(p/3)*(x*qbits*x+y*qbits*y+z*qbits*z);
}

