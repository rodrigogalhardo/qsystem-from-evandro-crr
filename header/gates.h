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
#include "using.h"
#include <map>
#include <armadillo>
#include <Python.h>

//! Class that holds some quantum gates
/*!
 * This class store the quantum gates of one qubit and the quantum gates
 * created by the user.
 */
class Gates {
  public:
    //! Constructor
    /*!
     * The constructor initialize the follow quantum gates:
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
     */
    Gates();

    //! Load constructor 
    /*!
     * Used to load gates from a file created by the method Gates::save.
     *
     * \param path to the file that store the quantum gates.
     *
     * \sa Gates::save
     */
    Gates(std::string path);

    //! Create an one qubit gate 
    /*!
     * The param `matrix` must have 4 elements organized like: `[a00, a01, a10,
     * a11]` = \f$\begin{bmatrix}a00&a01\\a10&a11\end{bmatrix}\f$.
     *
     * \param name of the new quantum gate.
     * \param matrix list of complex.
     *
     * \sa Gates::make_mgate Gates::make_cgate Gates::make_fgate
     */
    void make_gate(char name, vec_complex matrix);

    //! Create a quantum gate from a sparse matrix
    /*!
     * The matrix of the new quantum gate is created like: `U(row[i], col[i]) =
     * value[i]`.
     *
     * \param name of the new quantum gate.
     * \param size number of qubits affected by the new gate 
     * \param row list of row indices. 
     * \param col list of column indices 
     * \param value list of non-zero elements. 
     *
     * \sa Gates::make_gate Gates::make_cgate Gates::make_fgate
     */
    void make_mgate(std::string name,
                        size_t size,
                    vec_size_t row,
                    vec_size_t col,
                   vec_complex value);

    //! Create a controlled gate of X and Z
    /*!
     * Apply a sequence of gates `X`, `Z` and `I` if the `control` gates are in
     * the state \f$\left|1\right>\f$.
     *
     * \param name of the new quantum gate.
     * \param gates string with all the one qubit quantum gates, *e.g.* `"XZZXI"`.
     * \param control list of control gates.
     *
     * \sa Gates::make_gate Gates::make_mgate Gates::make_fgate
     */
    void make_cgate(std::string name,
                    std::string gates,
                     vec_size_t control);

    //! Create a quantum gate from a Python function
    /*! 
     * The Python function `func` must take as argumente and return an `int`.
     * This function must be defined in the range 0 to \f$2^\text{size}\f$-1.
     *
     * It's passible to pass an iterator that tells the order of creation of the matrix.
     *
     * \param name of the new quantum gate.
     * \param func Python function.
     * \param size number of qubits affected by the gate.
     * \param iterator with the matrix order of creation.
     *
     * \sa Gates::make_gate Gates::make_mgate Gates::make_cgate
     */
    void make_fgate(std::string name,
                      PyObject* func,
                         size_t size,
                      PyObject* iterator=Py_None);

    //! Get a string with information
    /*!
     *  This method is used in Python to cast a instance to `str`.
     *  
     *  \return String with information of the mutliple qubits gates.
     */
    std::string __str__();

    //! Save the multiple qubits quantum gates in a file
    /*!
     * The file created is a tar with all the multiple qubits quantum gates.
     *
     * \param path to the file that will be created.
     * \sa Gates::Gates
     */
    void save(std::string path);

    //! Return a quantum gate of one qubit
    /*!
     * This method is used by the QSystem class.
     * 
     * \param gate name of the gate.
     * \return Sparse matrix of the gate.
     * \sa Gates::mget
     */
    arma::sp_cx_mat& get(char gate);

    //! Retorn a quantum gate of multiple qubits
    /*!
     * This method is used by the QSystem class.
     * 
     * \param gate name of the gate.
     * \return Sparse matrix of the gate.
     * \sa Gates::get
     */
    arma::sp_cx_mat& mget(std::string gate);

  private:
  std::map<std::string, arma::sp_cx_mat> mmap;

  std::map<char, arma::sp_cx_mat> map{
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
};

