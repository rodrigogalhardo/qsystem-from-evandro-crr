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
    QSystem(size_t nqbits,
             Gates& gates,
            size_t seed=42,
       std::string state="vector");

    ~QSystem();
    
    /* evolution */
    void            evol(std::string gate,
                              size_t qbit, 
                              size_t count=1,
                                bool inver=false);
    void            cnot(size_t target, vec_size_t control);
    void            cphase(complex phase, size_t target, vec_size_t control);
    void            swap(size_t qbit_a, size_t qbit_b);
    void            qft(size_t qbegin, size_t qend, bool inver=false);
    
    /* measure */
    void             measure(size_t qbit, size_t count=1);
    void             measure_all();
    
    /* error channel */
    void            flip(char gate, size_t qbit, double p);
    void            amp_damping(size_t qbit, double p);
    void            dpl_channel(size_t qbit, double p);
    void            sum(size_t qbit, vec_str kraus, vec_float p);

    /* utility */
    std::string     __str__();
    size_t          size();
    std::string     state();

    void            save(std::string path);
    void            load(std::string path);

    void            change_to(std::string _state);

    vec_int         bits();
    PyObject*       get_qbits();
    void            set_qbits(vec_size_t row_ind,
                              vec_size_t col_ptr,
                             vec_complex values,
                                  size_t nqbits,
                             std::string state);

    /* ancilla */
    void            add_ancillas(size_t nqbits);
    void            rm_ancillas();

  private:
    void            sync();
    void            sync(size_t qbegin, size_t qend);
    void            clar();

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

    Gate_aux&       ops(size_t index);
    arma::sp_cx_mat get_gate(Gate_aux &op);
    cut_pair        cut(size_t &target, vec_size_t &control);
    void            fill(Gate_aux::Tag tag, size_t qbit, size_t size_n);
 

    Gates&           gates;
    size_t          _size;
    std::string     _state;
    Gate_aux*       _ops;
    bool            _sync;
    arma::sp_cx_mat qbits;
    Bit*            _bits;

    /* ancilla */
    size_t          an_size;
    Gate_aux*       an_ops;
    Bit*            an_bits;
};

