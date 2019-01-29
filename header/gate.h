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
#include <map>
#include <armadillo>
#include <vector>
#include <string>

using size_t = long unsigned;
using vec_cx = std::vector<std::complex<double>>;
using vec_size = std::vector<size_t>;

class Gate {
  public:
    Gate();

    Gate(std::string path);

    arma::sp_cx_mat& get(char gat);
    arma::sp_cx_mat& cget(std::string gat);

    void make_gate(char name,
                 vec_cx matrix);

    void make_gate(std::string name,
                        size_t size,
                      vec_size row,
                      vec_size col,
                        vec_cx value);

    void make_cgate(std::string name,
                    std::string gates,
                       vec_size control);

    void ls();

    void save(std::string path);

  private:
  std::map<std::string, arma::sp_cx_mat> cmap;

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

