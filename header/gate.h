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
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <string>

using size_t = long unsigned;
using cx = std::complex<double>;
using sp_cx_mat = Eigen::SparseMatrix<cx, Eigen::ColMajor>;
using sp_cx_vec = Eigen::SparseVector<cx>;
using cx_mat = Eigen::MatrixXcd;
using it_mat = sp_cx_mat::InnerIterator;
using it_vec = sp_cx_vec::InnerIterator;
using vec_cx = std::vector<cx>;
using vec_size = std::vector<size_t>;
/*
class mat_it {
    public:
      mat_it(sp_cx_mat& mat) : mat{mat} {}
      sp_cx_mat::InnerIterator& begin() {
        ncol = 1;
        curr = sp_cx_mat::InnerIterator(mat, 0);
        return curr;
      }
      sp_cx_mat::InnerIterator end() {
        if 
      }
      sp_cx_mat::InnerIterator& operator++() {
        ++curr;
        if (not curr and ncul < mat.outerSize()) {
          curr = sp_cx_mat::InnerIterator(mat, ncol)
          ++ncol;
        }
        return curr;
      }
    private:
      sp_cx_mat& mat;
      sp_cx_mat::InnerIterator curr;
      size_t ncol;  
};
*/ 
class Gate {
 public:
    Gate();

//  Gate(std::string path);

    sp_cx_mat& get(char gat);
    sp_cx_mat& cget(std::string gat);

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

    std::string __str__();

//  void save(std::string path);

  private:
  std::map<std::string, sp_cx_mat> cmap;

  std::map<char, sp_cx_mat> map;

};

