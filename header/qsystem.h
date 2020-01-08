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

#pragma once
#include "gate.h"
#include <variant>
#include <map>

//! Quantum circuit simulator class.
class QSystem {

  struct Gate_aux {
    Gate_aux();
    ~Gate_aux();

    bool busy();

    enum Tag {GATE_1, GATE_N, R,
              U3, CNOT, CPHASE,
              SWAP, QFT} tag;

    std::variant<char, r_pair, 
                 str, mat_ptr,
                 cnot_pair,
                 u3_tuple,
                 cph_tuple> data;

    size_t size;
    bool   invert;
  };

  enum Bit {NONE, ZERO, ONE};

  public:
    //! Constructor 
    /*!
     * The qubits are initialized in the state \f$\left|\text{init}\right>\f$.
     *
     * \param num_qbits number of qubits in the system.
     * \param seed for the pseudorandom number generator.
     * \param representation of the system, use `"bitwise"`, `"vector"` or
     * `"matrix"` (for density matrix).
     * \param init initial state.
     */
    QSystem(size_t num_qbits,
            size_t seed=42,
               str representation="bitwise",
            size_t init=0);

    ~QSystem();
    
    //! Apply a quantum gate
    /*!
     * The gate will be applied from the qubits `qbit` to `qbit+cout`. The
     * follow gates are available:
     * * ```'I'``` = \f$\begin{bmatrix}
     *                   1 & 0 \\
     *                   0 & 1
     *                  \end{bmatrix}\f$
     * * ```'Y'``` = \f$\begin{bmatrix}
     *                   0 & 1 \\
     *                   1 & 0
     *                  \end{bmatrix}\f$
     * * ```'X'``` = \f$\begin{bmatrix}
     *                   0 & -i \\
     *                   i & 0
     *                  \end{bmatrix}\f$
     * * ```'Z'``` = \f$\begin{bmatrix}
     *                   1 & 0 \\
     *                   0 & -1
     *                  \end{bmatrix}\f$
     * * ```'H'``` = \f${1\over\sqrt{2}}\begin{bmatrix}
     *                                   1 & 1 \\
     *                                   1 & -1
     *                                  \end{bmatrix}\f$
     * * ```'S'``` = \f$\begin{bmatrix}
     *                   1 & 0 \\
     *                   0 & i
     *                  \end{bmatrix}\f$
     * * ```'T'``` = \f$\begin{bmatrix}
     *                   1 & 0 \\
     *                   0 & e^{i\pi\over4}
     *                  \end{bmatrix}\f$.
     * 
     * \param gate name of the gate that will be user.
     * \param qbit qubit affected by the gate.
     * \param count number of successive repetitions of the gate.
     * \param invert if true, apply the inverse quantum gate.
     * \sa QSystem::rot QSystem::u3 QSystem::u2 QSystem::u1
     * QSystem::apply QSystem::cnot QSystem::cphase QSystem::qft QSystem::swap
     */
    void evol(char gate,
            size_t qbit, 
            size_t count=1,
              bool invert=false);

    //! Rotate a qubit in X, Y or Z axis
    /*! 
     * The rotation will be applied from the qubits `qbit` to `qbit+cout`
     * The rotation matrix for each axis are:
     * * ```'X'``` = \f$\begin{bmatrix}
     *                   \cos{\theta\over2} & -i\sin{\theta\over2} \\
     *                   -i\sin{\theta\over2} & \cos{\theta\over2} \\
     *                  \end{bmatrix}\f$
     * * ```'Y'``` = \f$\begin{bmatrix}
     *                   \cos{\theta\over2} & -\sin{\theta\over2} \\
     *                   \sin{\theta\over2} & \cos{\theta\over2} \\
     *                  \end{bmatrix}\f$
     * * ```'Z'``` = \f$\begin{bmatrix}
     *                   -e^{i{\theta\over2}} & 0 \\
     *                   0 & e^{i{\theta\over2}} \\
     *                  \end{bmatrix}\f$
     *
     * \param axis of the rotation
     * \param angle \f$\theta\f$ of the rotation
     * \param qbit qubit affected by the rotation.
     * \param count number of successive repetitions of the gate.
     * \sa QSystem::evol QSystem::u3 QSystem::u2 QSystem::u1
     * QSystem::apply QSystem::cnot QSystem::cphase QSystem::qft QSystem::swap
     */
    void rot(char axis,
           double angle,
           size_t qbit, 
           size_t count=1);

    //! Apply an arbitrary u3 gate
    /*!
     * The gate will be applied from the qubits `qbit` to `qbit+cout`.
     * \f[ 
     * u3 = 
     * \begin{bmatrix}
     * \cos{\theta\over2}          & -e^{i\lambda}\sin{\theta\over2} \\
     * e^{i\phi}\sin{\theta\over2} & e^{i(\lambda+\phi)}\cos{\theta\over2}
     * \end{bmatrix}
     * \f]
     *
     * \param theta = \f$\theta\f$
     * \param phi = \f$\phi\f$
     * \param lambd = \f$\lambda\f$
     * \param qbit qubit affected by the gate.
     * \param count number of successive repetitions of the gate.
     * \sa QSystem::evol QSystem::rot QSystem::u2 QSystem::u1
     * QSystem::apply QSystem::cnot QSystem::cphase QSystem::qft QSystem::swap
     */
    void u3(double theta,
            double phi,
            double lambd,
            size_t qbit, 
            size_t count=1);

    //! Apply an arbitrary u2 gate
    /*!
     * The gate will be applied from the qubits `qbit` to `qbit+cout`.
     * \f[ 
     * u2 = {1\over\sqrt{2}}
     * \begin{bmatrix}
     * 1 & -e^{i\lambda} \\
     * e^{i\phi} & e^{i(\lambda+\phi)}
     * \end{bmatrix}
     * \f]
     *
     * \param phi = \f$\phi\f$
     * \param lambd = \f$\lambda\f$
     * \param qbit qubit affected by the gate.
     * \param count number of successive repetitions of the gate.  
     * \sa QSystem::evol QSystem::rot QSystem::u3 QSystem::u1
     * QSystem::apply QSystem::cnot QSystem::cphase QSystem::qft QSystem::swap
     */
    void u2(double phi,
            double lambd, 
            size_t qbit, 
            size_t count=1);

    //! Apply an arbitrary u1 (phase) gate
    /*!
     * The gate will be applied from the qubits `qbit` to `qbit+cout`.
     * \f[ 
     * u1 = 
     * \begin{bmatrix}
     * 1 & 0 \\
     * 0 & e^{i\lambda}
     * \end{bmatrix}
     * \f]
     *
     * \param lambd = \f$\lambda\f$
     * \param qbit qubit affected by the gate.
     * \param count number of successive repetitions of the gate. 
     * \sa QSystem::evol QSystem::rot QSystem::u3 QSystem::u2
     * QSystem::apply QSystem::cnot QSystem::cphase QSystem::qft QSystem::swap
     */
    void u1(double lambd, 
            size_t qbit,
            size_t count=1);

    //! Apply a gate from a Gate class
    /*!
     * The gate will be applied from the qubits `qbit` to `qbit+cout*(size of
     * the gate)`.
     *
     * \param gate instace of Gate class
     * \param qbit qubit affected by the gate.
     * \param count number of successive repetitions of the gate. 
     * \param invert if true, apply the inverse quantum gate.
     * \sa QSystem::evol QSystem::rot QSystem::u3 QSystem::u2 QSystem::u1
     * QSystem::apply QSystem::cnot QSystem::cphase QSystem::qft QSystem::swap
     */
    void apply(Gate gate, size_t qbit, size_t count=1, bool invert=false);

    //! Apply a controlled not 
    /*!
     * Apply a not in the `target` qubit if all the `control` qubits are in the
     * representation \f$\left|1\right>\f$.
     *
     * \param target target qubit.
     * \param control list of control qubits.
     * \sa QSystem::evol QSystem::rot QSystem::u3 QSystem::u2 QSystem::u1
     * QSystem::cphase QSystem::qft QSystem::swap
     */
    void cnot(size_t target, vec_size_t control);

    //! Apply a controlled phase
    /*!
     * Apply \f$\begin{bmatrix}1&0\\0&e^\phi\end{bmatrix}\f$, where
     * \f$e^\phi\f$ = `phase`, in the `target` qubit if all the `control`
     * qubits are in the representation \f$\left|1\right>\f$.
     *
     * \param phase \f$e^\phi\f$ value.
     * \param target target qubit.
     * \param control list of control qubits.
     * \sa QSystem::evol QSystem::rot QSystem::u3 QSystem::u2 QSystem::u1
     * QSystem::apply QSystem::cnot QSystem::qft QSystem::swap
     */
    void cphase(complex phase, size_t target, vec_size_t control);

    //! Apply a quantum Fourier transformation
    /*!
     * Apply the QFT in the range of qubits (`qbit_begin`, `qbit_end`].
     *
     * \param qbit_begin first qubit affected. 
     * \param qbit_end last qubit affected +1.
     * \param invert if true, apply the inverse quantum gate.
     * \sa QSystem::evol QSystem::rot QSystem::u3 QSystem::u2 QSystem::u1
     * QSystem::apply QSystem::cnot QSystem::cphase QSystem::swap
     */
    void qft(size_t qbit_begin, size_t qbit_end, bool invert=false);

    //! Swap two qubit
    /*!
     * \param qbit_a qubit that gonna be swapped with `qbit_b`.
     * \param qbit_b qubit that gonna be swapped with `qbit_a`.
     * \sa QSystem::evol QSystem::rot QSystem::u3 QSystem::u2 QSystem::u1
     * QSystem::apply QSystem::cnot QSystem::cphase QSystem::qft 
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
     * The measurements results are stored in a list. Where the n-th item is
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
     * where \f$\sigma\f$ is \f$\begin{bmatrix}1&0\\0&1\end{bmatrix}\f$ if
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
     * that takes the qubit to the maximally mixed representation with probability `p`.
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

    //! Get system representation in a string
    /*!
     *  This method is used in Python to cast a instance to `str`.
     *
     * \return String with the system representation.
     * \sa QSystem::size QSystem::representation
     */
    str __str__();

    //! Get the number of qubits
    /*!
     * The ancillary qubits are include in the count.
     *
     * \return Number of qubits in the system.
     * \sa QSystem::representation 
     */
    size_t size();

    //! Get the system representation
    /*!
     * \return `"bitwise"`, `"vector"` or `"matrix"` (for density matrix).
     *
     * \sa QSystem::size QSystem::change_to
     */
    str representation();

    //! Save the quantum representation in a file
    /*!
     * The file is in a machine dependent binary format defined by the library
     * Armadillo.
     *
     * \param path to the file that will be created.
     * \sa QSystem::load
     */
    void save(str path);

    //! Load the quantum representation from a file
    /*!
     * When load, all qubits are set to non-ancillary.
     *
     * \param path to the file that will be loaded.
     * \sa QSystem::save
     */
    void load(str path);

    //! Change the system representation
    /*!
     * The change from vector representation to density matrix is done by
     * \f$\left|\psi\right> \rightarrow
     * \left|\psi\right>\mkern-7mu\left<\psi\right|\f$. It is not possible 
     * to change from density matrix representation to any other representation.
     * 
     * \param new_representation use "vector" to change vector representation a
     * \sa QSystem::representation
     */
    void change_to(str new_representation);

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
    py_obj get_qbits();

    //! Change the matrix of the quantum system
    /*!
     * This method is used in Python by the non-member function `set_matrix`.
     * ```python
     * def set_matrix(q, m):
     *     from scipy import sparse
     *     from math import log2
     *     m = sparse.csc_matrix(m)
     *     if m.shape[0] == m.shape[1]:
     *         representation = 'matrix'
     *     else:
     *         representation = 'vector'
     *     size = int(log2(m.shape[0]))
     *     q.set_qbits(m.indices.tolist(), m.indptr.tolist(), m.data.tolist(), size, representation)
     * ```
     * 
     * \param row_ind row indices.
     * \param col_ptr column pointers.
     * \param values non-zero values.
     * \param num_qbits number of qubits.
     * \param representation representation.
     * \sa QSystem::get_qbits
     */
    void set_qbits(vec_size_t row_ind,
                   vec_size_t col_ptr,
                  vec_complex values,
                       size_t num_qbits,
                          str representation);

    //! Add ancillary qubits
    /*!
     * The ancillaries qubits are added to the end of the system and can be used in 
     * any method.
     *
     * \param num_qbits number of ancillas added.
     * \param init initial state of the ancillas.
     * \sa QSystem::rm_ancillas
     */
    void add_ancillas(size_t num_qbits, size_t init=0);

    //! Remove all ancillary qubits
    /*!
     * If the state is in vector or bitwise representation the ancillas are measured
     * before been removed. If the state is in density matrix representation,
     * the ancillas are removed by a partial trace operation, without been
     * measured.
     *
     * \sa QSystem::rm_ancillas
     */
    void rm_ancillas();

    /* src/qs_evol.cpp */
    void            sync();
  private:
    void            sync(size_t qbit_begin, size_t qbit_end);
    Gate_aux&       ops(size_t index);
    arma::sp_cx_mat get_gate(Gate_aux &op);
    cut_pair        cut(size_t &target, vec_size_t &control);
    void            fill(Gate_aux::Tag tag, size_t qbit, size_t size_n);

    void            evol_x(size_t qbit);
    void            evol_y(size_t qbit);
    void            evol_z(size_t qbit);
    void            evol_s(size_t qbit, bool invert);
    void            evol_t(size_t qbit, bool invert);
    void            evol_h(size_t qbit);

    /* src/qs_make.cpp */
    arma::sp_cx_mat make_gate(arma::sp_cx_mat gate, size_t qbit);
    arma::sp_cx_mat make_rot(char axis, double angle);
    arma::sp_cx_mat make_cnot(size_t target,
                          vec_size_t control,
                              size_t size_n);
    arma::sp_cx_mat make_cphase(complex phase,
                                 size_t target,
                             vec_size_t control,
                                 size_t size_n);
    arma::sp_cx_mat make_swap(size_t size_n);
    arma::sp_cx_mat make_qft(size_t size_n);
    arma::sp_cx_mat make_u3(double theta, double phi, double lambd);

    /* src/qs_utility.cpp */
    void            clear();

    /*--------------------*/
    size_t          _size;
    str             _representation;
    Gate_aux*       _ops;
    bool            _sync;
    arma::sp_cx_mat qbits;
    dict            bwqbits;
    Bit*            _bits;

    size_t          an_size;
    Gate_aux*       an_ops;
    Bit*            an_bits;

    std::map<char, arma::sp_cx_mat> gates{
      {'I', arma::sp_cx_mat{arma::cx_mat{{{{1,0}, {0,0}},
                                          {{0,0}, {1,0}}}}}},
      {'X', arma::sp_cx_mat{arma::cx_mat{{{{0,0}, {1,0}},
                                          {{1,0}, {0,0}}}}}},
      {'Y', arma::sp_cx_mat{arma::cx_mat{{{{0,0}, {0,-1}},
                                          {{0,1}, {0,0}}}}}},
      {'Z', arma::sp_cx_mat{arma::cx_mat{{{{1,0}, {0,0}},
                                          {{0,0}, {-1,0}}}}}},
      {'H', arma::sp_cx_mat{(1/sqrt(2))*arma::cx_mat{{{{1,0}, {1,0}},
                                                      {{1,0}, {-1,0}}}}}},
      {'S', arma::sp_cx_mat{arma::cx_mat{{{{1,0}, {0,0}},
                                          {{0,0}, {0,1}}}}}},
      {'T', arma::sp_cx_mat{arma::cx_mat{{{{1,0}, {0,0}},
                                          {{0,0}, {1/sqrt(2),1/sqrt(2)}}}}}},
    };


    inline void valid_qbit(str name, size_t qbit);
    inline void valid_count(size_t qbit, size_t count, size_t size_n=1);
    inline void valid_control(vec_size_t &control);
    inline void valid_phase(complex phase);
    inline void valid_swap(size_t qbit_a, size_t qbit_b);
    inline void valid_range(size_t qbit_begin, size_t qbit_end);
    inline void valid_gate(str name, char gate);
    inline void valid_p(double p);
    inline void valid_state();
    inline void valid_krau(vec_str &kraus);
    inline void valid_state_str(str &representation);
    inline void valid_init(size_t init, size_t num_qbits);
    inline void valid_not_bw();

};

/******************************************************/
inline void QSystem::valid_qbit(str name, size_t qbit) {
  if (qbit >= size()) {
      sstr err;
      err << "\'" << name << "\' argument should be in the range of 0 to "
          << (size()-1);
      throw std::invalid_argument{err.str()};
  }
}

/******************************************************/
inline void QSystem::valid_count(size_t qbit, size_t count, size_t size_n) {
  if (count == 0 or qbit+count*size_n > size()) {
      sstr err;
      err << "\'cout\' argument should be greater than 0 "
          << "and \'qbit+count\' should be in the range of 0 to "
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
inline void QSystem::valid_range(size_t qbit_begin, size_t qbit_end) {
  if (qbit_begin >= size() or qbit_end > size() or qbit_begin >= qbit_end) {
      sstr err;
      err << "\'qbit_begin\' argument should be in the "
          << "range of 0 to " << (size()-1)
          << " and argument \'qbit_end\' should be greater than \'qbit_begin\' "
          << "and in the range of 1 to " << size();
      throw std::invalid_argument{err.str()};
  }
}

/******************************************************/
inline void QSystem::valid_gate(str name, char gate) {
  if (not (gate == 'X' or gate == 'Y' or gate == 'Z')) {
    sstr err;
    err << "\'" << name << "\' argument must be equal to \'X\', \'Y\' or \'X\'";
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
  if (_representation != "matrix") {
    sstr err;
    err << "\'representation\' must be in \"matrix\" to apply this channel";
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

inline void QSystem::valid_state_str(str &representation) {
  if (representation != "matrix" and representation != "vector" and representation != "bitwise") {
    sstr err;
    err << "\'representation\' argument must have value " 
        <<  "\"vector\", \"matrix\" or \"bitwise\", not \""
        << representation << "\"";
    throw std::invalid_argument{err.str()};
  }
}

inline void QSystem::valid_init(size_t init, size_t num_qbits) {
  if (init >= (1ul << num_qbits)) {
    sstr err;
    err << "\'init\' argument must be in the range of 0 to "
        << (1ul << num_qbits);
    throw std::invalid_argument{err.str()};
  }
}

inline void QSystem::valid_not_bw() {
  if (representation() == "bitwise") {
    sstr err;
    err << "\'representation\' can not be in \"bitwise\" to use this method";
    throw std::runtime_error{err.str()};
  }
}

