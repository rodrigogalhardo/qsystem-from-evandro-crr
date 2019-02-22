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
#include "gate.h"
#include <Python.h>
#include <variant>

class QSystem {

  struct Op {
    Op();
    ~Op();
    enum Tag {NONE, GATE_1,
              GATE_N, CNOT,
              CPHASE, SWAP, 
              QFT} tag;
    std::variant<char,
                 std::string,
                 cnot_pair,
                 cph_tuple> data;
    size_t size;
  };

  enum Bit {none, zero, one};

  public:
    QSystem(size_t nqbits,
            size_t seed,
             Gate& gate,
       std::string state);

    QSystem(std::string path,
                 size_t seed,
                  Gate& gate);
    ~QSystem();
    
    /* evolution */
    void            evol(char gate, size_t qbit);
    void            evol(char gate, size_t qbegin, size_t qend);
    void            evol(std::string gates);
    void            evol(std::string u, size_t qbit);
    void            cnot(size_t target, vec_size control);
    void            cphase(cx phase, size_t target, vec_size control);
    void            swap(size_t qbit_a, size_t qbit_b);
    void            qft(size_t qbegin, size_t qend);
    
    /* measure */
    void             measure(size_t qbit);
    void             measure(size_t qbegin, size_t qend);
    void             measure_all();
    
    /* error channel */
    void            flip(char gate, size_t qbit, double p);
    void            amp_damping(size_t qbit, double p);
    void            dpl_channel(size_t qbit, double p);
    void            sum(size_t qbit, vec_str kreaus, vec_d p);

    /* utility */
    std::string     __str__();
    size_t          get_size();
    std::string     get_state();
    void            save(std::string path);
    void            change_to(std::string state);
    vec_int         get_bits();
    PyObject*       get_qbits();
    void            set_qbits(vec_size row_ind,
                              vec_size col_ptr,
                                vec_cx values,
                                size_t nqbits,
                           std::string state);

    /* ancilla */
    void            add_ancillas(size_t an_num);
    void            rm_ancillas();
    void            an_evol(char gate, size_t qbit);
    void            an_evol(char gate, size_t qbegin, size_t qend);
    void            an_measure(size_t qbit);
    void            an_measure(size_t qbegin, size_t qend);
    size_t          get_an_size();
    vec_int         get_an_bits();

  private:
    void            sync();
    void            sync(size_t qbegin, size_t qend);
    void            clar();
    arma::sp_cx_mat make_gate(arma::sp_cx_mat gate, size_t qbit);
    arma::sp_cx_mat make_cnot(size_t target,
                            vec_size control,
                              size_t size_n);
    arma::sp_cx_mat make_cphase(cx phase,
                            size_t target,
                          vec_size control,
                            size_t size_n);
    arma::sp_cx_mat make_swap(size_t size_n);
    arma::sp_cx_mat make_qft(size_t size_n);
    arma::sp_cx_mat get_gate(Op &op);
    cut_pair        cut(size_t &target, vec_size &control);
    void            fill(Op::Tag tag, size_t qbit, size_t size_n);
 

    Gate&           gate;
    size_t          size;
    std::string     state;
    Op*             ops;
    bool            syncc;
    arma::sp_cx_mat qbits;
    Bit*            bits;

    /* ancilla */
    size_t          an_size;
    Op*             an_ops;
    Bit*            an_bits;
};

