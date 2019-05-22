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

namespace gate {
  class Gate {
    public:
      Gate(mat_ptr mat);
      Gate(str path);

      //! Get a string with information
      /*!
       *  This method is used in Python to cast a instance to `str`.
       *  
       *  \return String with information of the mutliple qubits gates.
       */
      str __str__();

      //! Save the multiple qubits quantum gates in a file
      /*!
       * The file created is a tar with all the multiple qubits quantum gates.
       *
       * \param path to the file that will be created.
       * \sa Gates::Gates
       */
      void save(str path);

      mat_ptr& get_mat();

    private:
      mat_ptr mat;
  };

  /*!
   * The param `matrix` must have 4 elements organized like: `[a00, a01, a10,
   * a11]` = \f$\begin{bmatrix}a00&a01\\a10&a11\end{bmatrix}\f$.
   *
   * \param matrix list of complex.
   *
   * \sa Gates::make_mgate Gates::make_cgate Gates::make_fgate
   */
  Gate make_gate(vec_complex matrix);

  //! Create a quantum gate from a sparse matrix
  /*!
   * The matrix of the new quantum gate is created like: `U(row[i], col[i]) =
   * value[i]`.
   *
   * \param size number of qubits affected by the new gate 
   * \param row list of row indices. 
   * \param col list of column indices 
   * \param value list of non-zero elements. 
   *
   * \sa Gates::make_gate Gates::make_cgate Gates::make_fgate
   */
  Gate make_mgate(size_t size,
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
   *
   * \sa Gates::make_gate Gates::make_mgate Gates::make_fgate
   */
  Gate make_cgate(str gates,
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
   *
   * \sa Gates::make_gate Gates::make_mgate Gates::make_cgate
   */
  Gate make_fgate(py_obj func,
                  size_t size,
                  py_obj iterator=Py_None);

}

