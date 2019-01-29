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
typedef long unsigned int size_t;
%}

%template(vec_size) std::vector<size_t>;
%template(vec_cx) std::vector<std::complex<double>>;
%template(vec_int) std::vector<int>;

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

class Code5:
    def __init__(self, gate, path='', size=0, seed=42, state='mix'):
        self.gate = gate;
        if path != '':
            self.q = QSystem(path, seed, self.gate)
            self.size = self.q.get_size()//5
        else:
            self.size = size;
            self.q = QSystem(5*size, seed, self.gate, state)
            self.prepare()

        self.syn = { 1:('Z',2),  2:('X',0),  3:('Z',3),  4:('X',3),
                     5:('X',1),  6:('Z',4),  7:('Y',3),  8:('Z',1),
                     9:('X',4), 10:('X',2), 11:('Y',2), 12:('Z',0),
                     13:('Y',1), 14:('Y',0), 15:('Y',4), 0:('I',0)}

        self.gate.make_gate('f', [-1,  0,
                                   0, -1])

    @staticmethod
    def make_gate(size):
        gate = Gate()
        pre_gate = ['XIZIZ', 'XZIZI', 'ZIZIX', 'IZIXZ']
        for pr in zip(pre_gate, [1, 4, 3, 2]):
            gate.make_cgate(pr[0], pr[0], [pr[1]])

        row, col, value = [], [], []
        for i in range(2**5):
            q0, q1 = int(bool(i & (1 << 4))), int(bool(i & (1 << 3)))
            q2, q3 = int(bool(i & (1 << 2))), int(bool(i & (1 << 1)))
            q4 = int(bool(i & 1))
            row.append((q1 << 4) | (q4 << 3) | (q2 << 2) | (q0 << 1) | q3)
            col.append(i)
            value.append(1)
        gate.make_gate('PERM_H', 5, row, col, value)
 
        ks = ['IZXXZ', 'ZIZXX', 'XZIZX', 'XXZIZ']
        for i in range(size):
            for k in ks:
                k_gate = k +  str().join(['IIIII' for _ in range(size-i-1)]) + 'I'
                gate.make_cgate(k + str(i), k_gate, [len(k_gate)-1])
        return gate

    def prepare(self):
        pre_gate = ['XIZIZ', 'XZIZI', 'ZIZIX', 'IZIXZ']

        for i in range(self.size):
            self.q.evol('Z', i*5)

        for pr in zip(pre_gate, [1, 4, 3, 2]):
            for i in range(self.size):
                self.q.evol('H', i*5+pr[1])
            for i in range(self.size):
                self.q.evol(pr[0], i*5)

    def correct(self, qbit):
        an_bits = []

        for k in ['IZXXZ', 'ZIZXX', 'XZIZX', 'XXZIZ']:
            self.q.add_ancillas(1)
            self.q.an_evol('H', 0)
            self.q.evol(k + str(qbit), qbit*5)
            self.q.an_evol('H', 0)
            self.q.an_measure(0)
            an_bits.append(self.q.get_an_bits()[0])
            self.q.rm_ancillas()
        i = sum([x[0]*x[1] for x in zip(an_bits, [1, 2, 4, 8])])
        if i != 0:
            error = self.syn[i]
            self.q.evol(error[0], qbit*5 + error[1])
        return an_bits, self.syn[i]

    def measure_all(self):
        self.q.measure_all()

    def measure(self, qbegin, qend):
        for i in range(qbegin, qend):
            self.q.measure(i*5, (i+1)*5)

    def H(self, qbit):
        self.q.evol('H', qbit*5, (qbit+1)*5)
        self.q.evol('PERM_H', qbit*5)
        self.q.evol('f', qbit*5, (qbit+1)*5)

    def print_state(self):
        self.q.print_state()

    def all_bits(self):
        return self.q.get_bits()

    def bit(self, i):
        return sum(self.q.get_bits()[i*5:(i+1)*5])%2

    def flip(self, sigma, lqbit, qbit, p):
        self.q.flip(sigma, lqbit*5+qbit, p)

    def evol(self, u, qbit):
        self.q.evol(u, 5*qbit, 5*(qbit+1))

    def cnot(self, target, control):
        for i in range(5):
            self.q.cnot(target*5+i, [control*5+i])

    def save(self, path):
        self.q.save(path)

%}

