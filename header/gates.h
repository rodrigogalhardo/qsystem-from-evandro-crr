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

class Gates {
  public:
    Gates();

    Gates(std::string path);

    arma::sp_cx_mat& get(char gate);
    arma::sp_cx_mat& mget(std::string gate);

    void make_gate(char name, vec_complex matrix);

    void make_mgate(std::string name,
                        size_t size,
                    vec_size_t row,
                    vec_size_t col,
                   vec_complex value);

    void make_cgate(std::string name,
                    std::string gates,
                     vec_size_t control);

    void make_fgate(std::string name,
                      PyObject* func,
                         size_t size,
                      PyObject* iterator=Py_None);

    std::string __str__();

    void save(std::string path);

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

