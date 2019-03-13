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

#include "../header/gates.h"
#include "../header/microtar.h"

using namespace arma;

/*********************************************************/
Gates::Gates() {}

/*********************************************************/
Gates::Gates(std::string path) {
  mtar_t tar;
  mtar_header_t h;
  char *m;

  mtar_open(&tar, path.c_str(), "r");

  for (; mtar_read_header(&tar, &h) != MTAR_ENULLRECORD; mtar_next(&tar)) {
    m = (char*) calloc(1, h.size);
    mtar_read_data(&tar, m, h.size);
    std::stringstream ss; 
    ss.write(m, h.size);
    free(m);
    sp_cx_mat matrix;
    matrix.load(ss, arma_binary);
    mmap[std::string(h.name)] = matrix;
  }

  mtar_close(&tar);
}

/*********************************************************/
sp_cx_mat& Gates::get(char gate) {
  return map.at(gate);
}

/*********************************************************/
sp_cx_mat& Gates::mget(std::string gate) {
  return mmap.at(gate);
}

/*********************************************************/
void Gates::make_gate(char name, vec_complex matrix) {
  map[name] = sp_cx_mat{cx_mat{{{matrix[0], matrix[1]},
                                {matrix[2], matrix[3]}}}};
}

/*********************************************************/
void Gates::make_mgate(std::string name,
                          size_t size, 
                      vec_size_t row,
                      vec_size_t col,
                     vec_complex value) {
  auto sizem = 1ul << size;
  sp_cx_mat m{sizem, sizem};

  for (size_t i = 0; i < row.size(); i++) {
    m(row[i], col[i]) = value[i];
  }

  mmap[name] = m;
}

/*********************************************************/
void Gates::make_cgate(std::string name,
                      std::string gates,
              std::vector<size_t> control) {

  size_t size = gates.size();
  size_t x = 0;
  size_t z = 0;
  for (size_t i = 0; i < size; i++) {
    if (gates[i] == 'X')
      x |= 1ul << (size-i-1);
    if (gates[i] == 'Z')
      z |= 1ul << (size-i-1);
  }

  auto parity = [&](size_t x) {
    for (int i = 32; i > 0; i /= 2)
      x ^= x >> i;
    return x & 1; 
  };

  sp_cx_mat cm{1ul << size, 1ul << size};

  for (size_t i = 0; i < (1ul << size); i++) {
    bool cond = true;
    for (size_t k = 0; k < control.size(); k++)
      cond = cond and ((1ul << (size-control[k]-1)) & i);

    if (cond) {
      size_t row = (i ^ x);
      cm(row, i) = pow(-1, parity(i & z));
    } else {
      cm(i,i) = 1;
    }
  }

  mmap[name] = cm;
}

/*********************************************************/
void Gates::make_fgate(std::string name,
                         PyObject* func,
                            size_t size,
                         PyObject* iterator) {

  sp_cx_mat m{1ul << size, 1ul << size};

  if (iterator == Py_None) {
    PyObject *builtins = PyEval_GetBuiltins(); 
    PyObject *range = PyDict_GetItemString(builtins , "range");
    iterator = PyEval_CallFunction(range, "i", 1ul << size);
  }

  auto* it = PyObject_GetIter(iterator);
  
  PyObject* pyj;
  while ((pyj = PyIter_Next(it))) {
    auto* arg = PyTuple_Pack(1, pyj);
    auto* pyi = PyObject_CallObject(func,   arg);
    
    auto i = PyLong_AsSize_t(pyi);
    auto j = PyLong_AsSize_t(pyj);

    m(i, j) = 1;
   
    Py_DECREF(pyi);
    Py_DECREF(pyj);
  }

  Py_DECREF(it);

  mmap[name] = m;
}

/*********************************************************/
std::string Gates::__str__() {
  std::stringstream out;
  for (auto& gate: mmap) {
    out << gate.first << " - "
        << log2(gate.second.n_rows)  << " qbits long"<< std::endl;
  }
  return out.str();
}

/*********************************************************/
void Gates::save(std::string path){
  mtar_t tar;
  mtar_open(&tar, path.c_str(), "w");

  for (auto &m : mmap) {
    std::stringstream file;
    m.second.save(file, arma_binary);
    file.seekg(0, ios::end);
    size_t size = file.tellg();
    file.seekg(0, ios::beg);
    mtar_write_file_header(&tar, m.first.c_str(), size);
    mtar_write_data(&tar, file.str().c_str(), size);
  }

  mtar_finalize(&tar);
  mtar_close(&tar);
}

