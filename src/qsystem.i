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

%include <std_string.i>
%include <std_vector.i>
%include <std_complex.i>

%inline %{
using size_t = long unsigned;
%}

%template(vec_size) std::vector<size_t>;
%template(vec_cx) std::vector<std::complex<double>>;
%template(vec_int) std::vector<int>;
%template(vec_d) std::vector<double>;
%template(vec_str) std::vector<std::string>;

%module qsystem
%{
  #include "../header/qsystem.h"
%}

%exception {
  try {
    $action
  } catch(std::exception &e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%typemap(out) std::vector<int> QSystem::get_bits %{
  $result = PyList_New($1.size());
  for (size_t i = 0; i < $1.size(); i++) {
    if ($1[i] == 0) 
      PyList_SetItem($result, i, Py_None);
    else 
      PyList_SetItem($result, i, PyLong_FromLong($1[i]-1));
  }
%}

%typemap(out) std::vector<int> QSystem::get_an_bits %{
  $result = PyList_New($1.size());
  for (size_t i = 0; i < $1.size(); i++) {
    if ($1[i] == 0) 
      PyList_SetItem($result, i, Py_None);
    else 
      PyList_SetItem($result, i, PyLong_FromLong($1[i]-1));
  }
%}

%include "../header/qsystem.h"
%include "../header/gate.h"

%pythoncode %{
def get_matrix(q):
    from scipy import sparse
    return sparse.csc_matrix(q.get_qbits()[0], q.get_qbits()[1])
 
def set_matrix(q, m, size, state):
    from scipy import sparse
    m = sparse.csc_matrix(m)
    q.set_qbits(m.indices.tolist(), m.indptr.tolist(), m.data.tolist(), size, state)

%}

