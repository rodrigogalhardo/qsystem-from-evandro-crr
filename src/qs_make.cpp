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
using namespace std::complex_literals;

/*********************************************************/
sp_cx_mat QSystem::make_gate(sp_cx_mat gate, size_t qbit) {
  sp_cx_mat m;
  size_t gate_size = log2(gate.n_rows);
  if (qbit == 0) {
    size_t eyesize = 1ul << (size()-gate_size);
    m = kron(gate, eye<sp_cx_mat>(eyesize, eyesize));
  } else if (qbit == size()-gate_size) {
    size_t eyesize = 1ul << (size()-gate_size);
    m = kron(eye<sp_cx_mat>(eyesize, eyesize), gate);
  } else {
    size_t eyesize = 1ul << qbit;
    m = kron(eye<sp_cx_mat>(eyesize, eyesize), gate);
    eyesize = 1ul << (size()-qbit-gate_size);
    m = kron(m, eye<sp_cx_mat>(eyesize, eyesize));
  }
  return m;
}

/*********************************************************/
sp_cx_mat QSystem::make_rot(char axis, double angle) {
  if (axis == 'X') 
    return sp_cx_mat{cx_mat{{{{cos(angle/2), 0}, {0, -sin(angle/2)}},
                             {{0, -sin(angle/2)}, {cos(angle/2), 0}}}}};
  else if (axis == 'Y')
    return sp_cx_mat{cx_mat{{{{cos(angle/2), 0},  {-sin(angle/2), 0}},
                             {{sin(angle/2), 0}, {cos(angle/2), 0}}}}};
  else 
    return sp_cx_mat{cx_mat{{{{-cos(angle/2), -sin(angle/2)}, {0, 0}},
                             {{0, 0}, {cos(angle/2), sin(angle/2)}}}}};
}

/*********************************************************/
sp_cx_mat QSystem::make_u3(double theta, double phi, double lambd) {
  return sp_cx_mat{cx_mat{{
    {{cos(theta/2), 0}, 
     {-cos(lambd)*sin(theta/2), -sin(lambd)*sin(theta/2)}},
    {{cos(phi)*sin(theta/2), sin(phi)*sin(theta/2)}, 
     {cos(lambd+phi)*cos(theta/2), sin(lambd+phi)*cos(theta/2)}
    }}}};
}

/*********************************************************/
sp_cx_mat QSystem::make_cnot(size_t target,
                         vec_size_t control,
                             size_t size_n) {
  sp_cx_mat cnotm{1ul << size_n, 1ul << size_n};

  for (size_t i = 0; i < (1lu << size_n); i++) {
    bool cond = true;

    for (size_t k = 0; (k < control.size()) and cond; k++)
      cond = cond and ((1ul << (size_n-control[k]-1)) & i);

    if (cond)
      cnotm(i, i ^ (1ul  << (size_n-target-1))) = 1; 
    else 
      cnotm(i, i) = 1;
  }

  return cnotm;
}

/******************************************************/
sp_cx_mat QSystem::make_cphase(complex phase,
                                size_t target,
                            vec_size_t control,
                                size_t size_n) {
  sp_cx_mat cphasem{1ul << size_n, 1ul << size_n};

  for (size_t i = 0; i < (1lu << size_n); i++) {
    bool cond = true;

    for (size_t k = 0; (k < control.size()) and cond; k++)
      cond = cond and ((1ul << (size_n-control[k]-1)) & i);

    cphasem(i, i) = cond and (i & (1ul << (size_n-target-1)))? phase : 1; 
  }

  return cphasem;
}

/******************************************************/
sp_cx_mat QSystem::make_swap(size_t size_n) {
  sp_cx_mat swapm{1ul << size_n, 1ul << size_n};

  for (size_t i = 0; i < (1ul << (size_n-1)); i++) {
    if (i%2 == 1) 
      swapm((i | (1ul << (size_n-1))) ^ 1ul, i) = 1;
    else 
      swapm(i, i) = 1;
  }

  for (size_t i = 0; i < (1ul << (size_n-1)); i++) {
    if (i%2 == 0) 
      swapm(i ^ 1ul, i | (1ul << (size_n-1))) = 1;
    else 
      swapm(i | (1ul << (size_n-1)), i | (1ul << (size_n-1))) = 1;

  }

  return swapm;
}

/******************************************************/
sp_cx_mat QSystem::make_qft(size_t size_n) {
  double pi = acos(-1);
  complex w = std::exp((2*pi*1i)/(pow(2, size_n)));
  sp_cx_mat qftm{1ul << size_n, 1ul << size_n};
  for (size_t i = 0; i < (1ul << size_n); i++) 
    for (size_t j = 0; j < (1ul << size_n); j++) 
      qftm(i, j) = (1/sqrt(1 << size_n))*pow(w, i*j);
  return qftm;
}

