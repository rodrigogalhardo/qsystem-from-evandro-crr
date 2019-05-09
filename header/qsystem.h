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

//! Quantum circuit simulator class.
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
    //! Constructor 
    /*!
     * All qubits are initialized in the state \f$\left|0\right>\f$.
     *
     * \param nqbits number of qubits in the system.
     * \param gates instance of class Gates that holds the gates used in the
     * method QSystem::evol.
     * \param seed for the pseudorandom number generator.
     * \param state representation of the system, use `"vector"` for vector.
     * state and `"matrix"` for density matrix
     */
    QSystem(size_t nqbits,
             Gates& gates,
            size_t seed=42,
       std::string state="vector");

    ~QSystem();
    
    //! Apply a quantum gate
    /*!
     * The gate will be applied in from qubits `qbit` to `qbit*cout*(size of
     * the gate)`. If `gate` parameter is just one character long, the size of
     * the gate is necessarily one.
     *
     * \param gate name of the gate that will be user.
     * \param qbit qubit affected by the gate.
     * \param count number of successive repetitions of the gate.
     * \param inver if true, apply the inverse quantum gate.
     * \sa QSystem::cnot QSystem::cphase QSystem::qft QSystem::swap
     */
    void evol(std::string gate,
                   size_t qbit, 
                   size_t count=1,
                     bool inver=false);
    //! Apply a controlled not 
    /*!
     * Apply a not in the `target` qubit if all the `control` qubits are in the
     * state \f$\left|1\right>\f$.
     *
     * \param target target qubit.
     * \param control list of control qubits.
     * \sa QSystem::evol QSystem::cphase QSystem::qft QSystem::swap
     */
    void cnot(size_t target, vec_size_t control);

    //! Apply a controlled phase
    /*!
     * Apply \f$\begin{bmatrix}1&0\\0&e^\phi\end{bmatrix}\f$, whare
     * \f$e^\phi\f$ = `phase`, in the `target` qubit if all the `control`
     * qubits are in the state \f$\left|1\right>\f$.
     *
     * \param phase \f$e^\phi\f$ value.
     * \param target target qubit.
     * \param control list of control qubits.
     * \sa QSystem::evol QSystem::cnot QSystem::qft QSystem::swap
     */
    void cphase(complex phase, size_t target, vec_size_t control);

    //! Apply a quantum Fourier transformation
    /*!
     * Apply the QFT in the range of qubits (`qbegin`, `qend`].
     *
     * \param qbegin first qubit affected. 
     * \param qend last qubit affected +1.
     * \sa QSystem::evol QSystem::cnot QSystem::cphase QSystem::swap
     */
    void qft(size_t qbegin, size_t qend, bool inver=false);

    //! Swap two qubit
    /*!
     * \param qbit_a qubit that gonna be swapped with `qbit_b`.
     * \param qbit_a qubit that gonna be swapped with `qbit_a`.
     * \sa QSystem::evol QSystem::cnot QSystem::cphase QSystem::qft
     */
    void swap(size_t qbit_a, size_t qbit_b);
    
    //! Measure qubits in the computational base
    /*!
     * Measure qubits from `qbit` to `qbit*cout`. All the measurements results
     * are assessable throw the QSystem::bits method.
     *
     * \param qbit qubit affected by the measurement.
     * \param count number qubits measured from `qbit`.
     * \sa QSystem::measure_all QSystem::bits
     */
    void measure(size_t qbit, size_t count=1);

    //! Measure all qubits in the computational base
    /*!
     * The measurements results are assessable throw the QSystem::bits method.
     * \sa QSystem::measure QSystem::bits
     */
    void measure_all();

    //! Get the measurements results
    /*! 
     * The measurements results are stored in a list. Whare the n-th item is
     * the measurement result of the qubit n. If the n-th qubits has never been
     * measured it's value is `None`.
     *
     * \return List of the measurement result.
     * \sa QSystem::measure QSystem::measure_all
     */
    vec_int bits();
    
    //! Apply a bit, phase or bit-phase flip error
    /*!
     * Apply the Kraus operator
     * \f[
     *    E_0 = \sqrt{p}\sigma\\
     *    E_1 = \sqrt{1-p}I,
     * \f]
     * whare \f$\sigma\f$ is \f$\begin{bmatrix}1&0\\0&1\end{bmatrix}\f$ if
     * `gate` is ```'X'```,  \f$\begin{bmatrix}1&0\\0&1\end{bmatrix}\f$ if
     * `gate` is ```'Z'``` or \f$\begin{bmatrix}1&0\\0&-1\end{bmatrix}\f$ if
     * `gate` is ```'Y'```.
     * 
     * \param gate use ```'X'``` for bit flip, ```'Z'``` for phase flip 
     * or ```'Y'``` for bit-phase flip.
     * \param qbit qubit effected by the error.
     * \param p probability of the error occur.
     * \sa QSystem::amp_damping QSystem::dpl_channel QSystem::sum
     */
    void flip(char gate, size_t qbit, double p);

    //! Apply an amplitude damping channel error
    /*!
     * Apply the Kraus operator
     * \f[
     *    E_0 = \begin{bmatrix}1&0\\0&\sqrt{1-p}\end{bmatrix}\\
     *    E_1 = \begin{bmatrix}0&0\\ \sqrt{p}&0\end{bmatrix},
     * \f]
     *
     * The system must be in density matrix representation to use this method,
     * otherwise you can use the flow code to achieve a similar result:
     * ```python
     * def amp_damping(q, qbit, p):
     *    from random import choices
     *    if choices([True, False], weights=[p, 1-p])[0]:
     *        q.measure(qbit)
     *        if q.bits[qbit] == 1:
     *            q.evol('X', qbit)
     * ```
     *
     * \sa QSystem::flip QSystem::dpl_channel QSystem::sum
     */
    void amp_damping(size_t qbit, double p);

    //! Apply a depolarization channel error
    /*!
     * Apply the operator
     * \f[
     *    \mathcal{E}(\rho) = \left(1-{3p\over4}\right)\rho
     *    +{p\over4}(X\rho X+Y\rho Y+Z\rho Z),
     * \f]
     * that takes the qubit to the maximally mixed state with probability `p`.
     *
     * The system must be in density matrix representation to use this method.
     *
     * \param qbit qubit effected by the error.
     * \param p probability of the error occur.
     * \sa QSystem::flip QSystem::amp_damping QSystem::sum
     */
    void dpl_channel(size_t qbit, double p);
  
    //! Apply a sum operator
    /*!
     * To apply some Kraus operator like
     * \f[
     *    E_1 = {\sqrt{p_1}} (U_{11}\otimes\dots\otimes U_{1n})\\
     *    E_2 = {\sqrt{p_2}} (U_{21}\otimes\dots\otimes U_{2n})\\
     *    \vdots\\
     *    E_m = {\sqrt{p_m}} (U_{m1}\otimes\dots\otimes U_{mn}),
     * \f]
     * pass the follow parameters
     * * `kraus` = [\f$U_{11}\otimes\dots\otimes U_{1n},\,
     * U_{21}\otimes\dots\otimes U_{2n},\, \dots,\, U_{m1}\otimes\dots\otimes
     * U_{mn}\f$] and 
     * * `p` = [\f$p_1,\,p_2,\,\dots,\,p_m\f$]
     *
     * The system must be in density matrix representation to use this method,
     * otherwise you can use the flow code to achieve a similar result:
     * ```python
     * def sum(q, qbit, kraus, p):
     *     from random import choices
     *     aux = 0
     *     for gate in choices(kraus, weights=p)[0]:
     *         q.evol(gate, qbit+aux)
     *         aux += 1
     * ```
     *
     * \param qbit first qubit effected by the error.
     * \param kraus Kraus operators list.
     * \param p probability list.
     * \sa QSystem::flip QSystem::amp_damping QSystem::dpl_channel
     */
    void sum(size_t qbit, vec_str kraus, vec_float p);

    //! Get system state in a string
    /*!
     *  This method is used in Python to cast a instance to `str`.
     *
     * \return String with the system state.
     * \sa QSystem::size QSystem::state
     */
    std::string __str__();

    //! Get the number of qubits
    /*!
     * The ancillary qubits are include in the count.
     *
     * \return Number of qubits in the system.
     * \sa QSystem::state 
     */
    size_t size();

    //! Get the system representation
    /*!
     * \return `"vector"` for vector representation and `"matrix"` for density
     * matrix.
     *
     * \sa QSystem::size QSystem::change_to
     */
    std::string state();

    //! Save the quantum state in a file
    /*!
     * The file is in a machine dependent binary format defined by the library
     * Armadillo.
     *
     * \param path to the file that will be created.
     * \sa QSystem::load
     */
    void save(std::string path);

    //! Load the quantum state from a file
    /*!
     * When load, all qubits are set to non-ancillary.
     *
     * \param path to the file that will be loaded.
     * \sa QSystem::save
     */
    void load(std::string path);

    //! Change the system representation
    /*!
     * The change from vector representation to density matrix is done by
     * \f$\left|\psi\right> \rightarrow
     * \left|\psi\right>\mkern-7mu\left<\psi\right|\f$. But, in the change from
     * density matrix to vector representation, just the measurement
     * probability is maintained. 
     *
     * \param new_state use "vector" to change vector representation a
     * \sa QSystem::state
     */
    void change_to(std::string new_state);

    //! Get the matrix of the quantum system 
    /*!
     * This method is used in Python by the non-member function `get_matrix`.
     * ```python
     * def get_matrix(q):
     *    from scipy import sparse
     *    return sparse.csc_matrix(q.get_qbits()[0], q.get_qbits()[1])
     * ```
     *
     * \return Tuple used to initialize a scipy space matrix.
     * \sa QSystem::set_qbits
     */
    PyObject* get_qbits();

    //! Change the matrix of the quantum system
    /*!
     * This method is used in Python by the non-member function `set_matrix`.
     * ```python
     * def set_matrix(q, m):
     *     from scipy import sparse
     *     from math import log2
     *     m = sparse.csc_matrix(m)
     *     if m.shape[0] == m.shape[1]:
     *         state = 'matrix'
     *     else:
     *         state = 'vector'
     *     size = int(log2(m.shape[0]))
     *     q.set_qbits(m.indices.tolist(), m.indptr.tolist(), m.data.tolist(), size, state)
     * ```
     * 
     * \param row_ind row indices.
     * \param col_ptr column pointers.
     * \param value non-zero values.
     * \param nqbits number of qubits.
     * \param state representation.
     * \sa QSystem::get_qbits
     */
    void set_qbits(vec_size_t row_ind,
                   vec_size_t col_ptr,
                  vec_complex values,
                       size_t nqbits,
                  std::string state);

    //! Add ancillary qubits
    /*!
     * The ancillaries qubits are added to the end of the system and can be used in 
     * any method.
     *
     * \param nqbits number of ancillas added.
     * \sa QSystem::rm_ancillas
     */
    void add_ancillas(size_t nqbits);

    //! Remove all ancillary qubits
    /*!
     * If the state is in vector representation the ancillas are measured
     * before been removed. If the state is in density matrix representation,
     * the ancillas are removed by a partial trace operation, without been
     * measured.
     *
     * \param nqbits number of ancillas added.
     * \sa QSystem::rm_ancillas
     */
    void rm_ancillas();

  private:
    /* src/qs_evol.cpp */
    void            sync();
    void            sync(size_t qbegin, size_t qend);
    Gate_aux&       ops(size_t index);
    arma::sp_cx_mat get_gate(Gate_aux &op);
    cut_pair        cut(size_t &target, vec_size_t &control);
    void            fill(Gate_aux::Tag tag, size_t qbit, size_t size_n);

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
    void            clear();

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
  if (qbegin >= size() or qend > size() or qbegin >= qend) {
      sstr err;
      err << "\'qbegin\' argument should be in the "
          << "range of 0 to " << (size()-1)
          << " and argument \'qend\' should be greater than \'qbegin\' "
          << "and in the range of 1 to " << size();
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

