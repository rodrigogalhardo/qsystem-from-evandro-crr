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
#include "using.h"

//! Quantum gate holder class
class Gate {
  public:
    //! Load constructor
    /*!
     * Load a quantum gate from a file.
     *
     * \param path to the file that will be loaded.
     */
    Gate(str path);

    //! Constructor 
    /*!
     * Constructor used by the functions from_matrix, from_sp_matrix,
     * cxz_gate and from_func.
     *
     * \param mat sparse matrix pointer.
     * \sa Gate::from_matrix Gate::from_sp_matrix Gate::cxz_gate Gate::from_func
     */
    Gate(mat_ptr mat, set_mat bwgate);

    //! Get a string with the matrix
    /*!
     *  This method is used in Python to cast a instance to `str`.
     *  
     *  \return String with the matrix.
     */
    str __str__();

    //! Save the quantum gate in a file
    /*!
     * The file is in a machine dependent binary format defined by the library
     * Armadillo.
     *
     * \param path to the file that will be created.
     * \sa Gate::Gate(str)
     */
    void save(str path);

    //! Create an one qubit gate 
    /*!
     * The param `matrix` must have 4 elements organized like: `[a00, a01, a10,
     * a11]` = \f$\begin{bmatrix}a00&a01\\a10&a11\end{bmatrix}\f$.
     *
     * \param matrix list of complex.
     * \return Gate class.
     * \sa Gate::from_sp_matrix Gate::cxz_gate Gate::from_func
     */
    static Gate from_matrix(vec_complex matrix);

    //! Create a quantum gate from a sparse matrix
    /*!
     * The matrix of the new quantum gate is created like: `U(row[i], col[i]) =
     * value[i]`.
     *
     * \param size number of qubits affected by the new gate 
     * \param row list of row indices. 
     * \param col list of column indices 
     * \param value list of non-zero elements. 
     * \return Gate class.
     * \sa Gate::from_matrix Gate::cxz_gate Gate::from_func
     */
    static Gate from_sp_matrix(size_t size,
                           vec_size_t row,
                           vec_size_t col,
                          vec_complex value);

    //! Create a controlled gate of X and Z
    /*!
     * Apply a sequence of gates `X`, `Z` and `I` if the `control` gates are in
     * the state \f$\left|1\right>\f$.
     *
     * \param gates string with all the one qubit quantum gates, *e.g.* `"XZZXI"`.
     * \param control list of control gates.
     * \return Gate class.
     * \sa Gate::from_matrix Gate::from_sp_matrix Gate::from_func
     */
    static Gate cxz_gate(str gates,
                  vec_size_t control);

    //! Create a quantum gate from a Python function
    /*! 
     * The Python function `func` must take as argumente and return an `int`.
     * This function must be defined in the range 0 to \f$2^\text{size}\f$-1.
     *
     * It's passible to pass an iterator that tells the order of creation of the matrix.
     *
     * \param func Python function.
     * \param size number of qubits affected by the gate.
     * \param iterator with the matrix order of creation.
     * \return Gate class.
     * \sa Gate::from_matrix Gate::from_sp_matrix Gate:cxz_gate
     */
    static Gate from_func(py_function func,
                               size_t size,
                          py_iterator iterator=Py_None);

    //! Return the matrix pointer
    /*! 
     * Method used by the QSystem class.
     * 
     * \return matrix pointer.
     * \sa QSystem
     */
    mat_ptr& get_mat();

    //! Return the bitwise gate pointer
    /*! 
     * Method used by the QSystem class.
     * 
     * \return gate pointer.
     * \sa QSystem
     */
    set& get_bwgate(size_t i);

    private:
    mat_ptr mat;
    set_mat bwgate;
};

