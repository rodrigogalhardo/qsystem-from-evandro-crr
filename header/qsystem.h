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

#pragma once
#include "gate.h"
#include <Python.h>

enum Bit {none, zero, one};
using vec_str = std::vector<std::string>;
using vec_int = std::vector<int>;

class QSystem {
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
    void            cnot(size_t target, vec_size control);
    void            evol(std::string u, size_t qbit);
    
    /* measure */
    void             measure(size_t qbit);
    void             measure(size_t qbegin, size_t qend);
    void             measure_all();
    
    /* error channel */
    void            flip(char gate, size_t qbit, double p);
    void            amp_damping(size_t qbit, double p);
    void            dpl_channel(size_t qbit, double p);

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
    void            clar();
    arma::sp_cx_mat make_gate(arma::sp_cx_mat gate, size_t qbit);

    Gate&           gate;
    size_t          size;
    std::string     state;
    char*           ops;
    vec_str         mops;
    bool            syncc;
    arma::sp_cx_mat qbits;
    Bit*            bits;

    /* ancilla */
    size_t          an_size;
    char*           an_ops;
    Bit*            an_bits;
};

