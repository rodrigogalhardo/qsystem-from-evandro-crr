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

#pragma once
#include "gates.h"
#include <Python.h>
#include <variant>

class QSystem {

  /* src/qs_utility.cpp */
  struct Gate_aux {
    Gate_aux();
    ~Gate_aux();

    bool busy();

    enum Tag {GATE_1, GATE_N,
              CNOT, CPHASE,
              SWAP, QFT} tag;

    std::variant<char,
                 std::string,
                 cnot_pair,
                 cph_tuple> data;

    size_t size;
    bool   inver;
  };

  enum Bit {NONE, ZERO, ONE};

  public:
    /* src/qs_utility.cpp */
    QSystem(size_t nqbits,
             Gates& gates,
            size_t seed=42,
       std::string state="vector");

    ~QSystem();
    
    /* src/qs_evol.cpp */
    void            evol(std::string gate,
                              size_t qbit, 
                              size_t count=1,
                                bool inver=false);
    void            cnot(size_t target, vec_size_t control);
    void            cphase(complex phase, size_t target, vec_size_t control);
    void            swap(size_t qbit_a, size_t qbit_b);
    void            qft(size_t qbegin, size_t qend, bool inver=false);
    
    /* src/qs_measure.cpp */
    void             measure(size_t qbit, size_t count=1);
    void             measure_all();
    vec_int          bits();
    
    /* src/qs_errors.cpp */
    void            flip(char gate, size_t qbit, double p);
    void            amp_damping(size_t qbit, double p);
    void            dpl_channel(size_t qbit, double p);
    void            sum(size_t qbit, vec_str kraus, vec_float p);

    /* src/qs_utility.cpp */
    std::string     __str__();
    size_t          size();
    std::string     state();

    void            save(std::string path);
    void            load(std::string path);

    void            change_to(std::string new_state);

    PyObject*       get_qbits();
    void            set_qbits(vec_size_t row_ind,
                              vec_size_t col_ptr,
                             vec_complex values,
                                  size_t nqbits,
                             std::string state);

    /* src/qs_ancillas.cpp */
    void            add_ancillas(size_t nqbits);
    void            rm_ancillas();

  private:
    /* src/qs_evol.cpp */
    void            sync();
    void            sync(size_t qbegin, size_t qend);
    void            clear();

    /* src/qs_make.cpp */
    arma::sp_cx_mat make_gate(arma::sp_cx_mat gate, size_t qbit);
    arma::sp_cx_mat make_cnot(size_t target,
                          vec_size_t control,
                              size_t size_n);
    arma::sp_cx_mat make_cphase(complex phase,
                                 size_t target,
                             vec_size_t control,
                                 size_t size_n);
    arma::sp_cx_mat make_swap(size_t size_n);
    arma::sp_cx_mat make_qft(size_t size_n);

    /* src/qs_utility.cpp */
    Gate_aux&       ops(size_t index);
    arma::sp_cx_mat get_gate(Gate_aux &op);
    cut_pair        cut(size_t &target, vec_size_t &control);
    void            fill(Gate_aux::Tag tag, size_t qbit, size_t size_n);

    /*--------------------*/
    Gates&           gates;
    size_t          _size;
    std::string     _state;
    Gate_aux*       _ops;
    bool            _sync;
    arma::sp_cx_mat qbits;
    Bit*            _bits;

    size_t          an_size;
    Gate_aux*       an_ops;
    Bit*            an_bits;

    /* src/qs_valid.cpp */
    inline void     valid_qbit(std::string name, size_t qbit);
    inline void     valid_count(size_t qbit, size_t count, size_t size_n=1);
    inline void     valid_control(vec_size_t &control);
    inline void     valid_phase(complex phase);
    inline void     valid_swap(size_t qbit_a, size_t qbit_b);
    inline void     valid_range(size_t qbegin, size_t qend);
    inline void     valid_gate(char gate);
    inline void     valid_p(double p);
    inline void     valid_state();
    inline void     valid_krau(vec_str &kraus);

};

/******************************************************/
inline void QSystem::valid_qbit(std::string name, size_t qbit) {
  if (qbit >= size()) {
      sstr err;
      err << "\'" << name << "\' argument should be in the range of 0 to "
          << (size()-1);
      throw std::invalid_argument{err.str()};
  }
}

/******************************************************/
inline void QSystem::valid_count(size_t qbit, size_t count, size_t size_n) {
  if (count == 0 and qbit+count*size_n <= size()) {
      sstr err;
      err << "\'cout\' argument should be greater than 0 "
          << "and \'qbit+count\' suld be in the range of 0 to "
          << size();
      throw std::invalid_argument{err.str()};
  }
}

/******************************************************/
inline void QSystem::valid_control(vec_size_t &control) {
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
inline void QSystem::valid_phase(complex phase) {
  if (std::abs(std::abs(phase) - 1.0) > 1e-14) {
    sstr err;
    err << "abs(phase) must be equal to 1";
    throw std::invalid_argument{err.str()};
  }
}

/******************************************************/
inline void QSystem::valid_swap(size_t qbit_a, size_t qbit_b) {
  if (qbit_a >= size() or qbit_b >= size()) {
      sstr err;
      err << "Arguments \'qbit_a\' and \'qbit_b\' should be in the "
          << "range of 0 to " << (size()-1);
      throw std::invalid_argument{err.str()};
  }
}

/******************************************************/
inline void QSystem::valid_range(size_t qbegin, size_t qend) {
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
inline void QSystem::valid_gate(char gate) {
  if (not (gate == 'X' or gate == 'Y' or gate == 'Z')) {
    sstr err;
    err << "\'gate\' argument must be equal to \'X\', \'Y\' or \'X\'";
    throw std::invalid_argument{err.str()};    
  }
}

/******************************************************/
inline void QSystem::valid_p(double p) {
  if (p < 0 or p > 1) {
    sstr err;
    err << "\'p\' argument should be in the range of 0.0 to 1.0";
    throw std::invalid_argument{err.str()};
  }
}

inline void QSystem::valid_state() {
  if (_state == "vector") {
    sstr err;
    err << "\'state\' must be in \"matrix\" to apply this channel";
    throw std::runtime_error{err.str()};
  }
}

inline void QSystem::valid_krau(vec_str &kraus) {
  size_t ksize = kraus[0].size();
  for (auto& k : kraus) {
    if (k.size() != ksize) {
      sstr err;
      err << "All \'kraus\' operators must have the same size";
      throw std::runtime_error{err.str()};
    }
  }
}

