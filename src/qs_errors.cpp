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
void QSystem::flip(char gate, size_t qbit, double p) {
  valid_gate("gate", gate);
  valid_qbit("qbit", qbit);
  valid_p(p);

  if (_representation == "vector") {
    if (auto pr = double(std::rand())/double(RAND_MAX); p != 0 and pr <= p) 
      evol(gate, qbit);

  } else if (_representation == "matrix") {
    sync();
    sp_cx_mat E0 = make_gate(gates[gate], qbit)*sqrt(p);

    size_t eyesize = 1ul << size();
    sp_cx_mat E1 = eye<sp_cx_mat>(eyesize, eyesize)*sqrt(1.f-p);

    qbits = E0*qbits*E0.t() + E1*qbits*E1;
  }
}

/******************************************************/
void QSystem::amp_damping(size_t qbit, double p) {
  valid_state();
  valid_qbit("qbit", qbit);
  valid_p(p);

  sync();

  sp_cx_mat E0 = make_gate(sp_cx_mat{cx_mat{{{{1, 0}, {0, 0}},
                                             {{0, 0}, {sqrt(1-p), 0}}}}}, qbit);
  sp_cx_mat E1 = make_gate(sp_cx_mat{cx_mat{{{{0, 0}, {sqrt(p), 0}},
                                             {{0, 0}, {0, 0}}}}}, qbit);
  qbits = E0*qbits*E0 + E1*qbits*E1.t();
}

/******************************************************/
void QSystem::dpl_channel(size_t qbit, double p) {
  valid_state();
  valid_qbit("qbit", qbit);
  valid_p(p);

  sync();

  sp_cx_mat X = make_gate(gates['X'], qbit);
  sp_cx_mat Y = make_gate(gates['Y'], qbit);
  sp_cx_mat Z = make_gate(gates['Z'], qbit);
  
  qbits = (1-p)*qbits+(p/3)*(X*qbits*X+Y*qbits*Y.t()+Z*qbits*Z);
}

/******************************************************/
void QSystem::sum(size_t qbit, vec_str kraus, vec_float p) {
  valid_state();
  valid_qbit("qbit", qbit);
  valid_krau(kraus);
    
  sync();

  sp_cx_mat qbits_tmp{qbits.n_rows, qbits.n_cols};

  for (size_t i = 0; i < kraus.size(); ++i) {
    sp_cx_mat E = gates[kraus[i][0]];
    for (size_t j = 1; j < kraus[i].size(); ++j)
      E = kron(E, gates[kraus[i][j]]);
    E = make_gate(E, qbit);
    
    qbits_tmp +=  p[i]*(E*qbits*E.t());
  }

  qbits = qbits_tmp;  
}

