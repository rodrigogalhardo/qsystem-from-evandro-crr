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

#include "../header/qsystem.h"
#include <iomanip>
#include <sstream>

using namespace arma;

QSystem::QSystem(size_t nqbits, size_t seed, Gate& gate, std::string state) :
  gate{gate}, size{nqbits}, state{state}, ops{new char[size]}, syncc{true},
  qbits{1lu << size, state == "mix" ? 1lu << size : 1},
  bits{new Bit[size]()}, an_size{0}, an_ops{nullptr}, an_bits{nullptr}
{
  if (state != "mix" and state != "pure") 
    throw std::invalid_argument{"Argument \'state\' must be \"pure\" or \"mix\", not \""
      + state + "\"."};
  qbits(0,0) = 1;
  std::memset(ops, 'I', size*sizeof(char));
  std::srand(seed);
}

QSystem::QSystem(std::string path, size_t seed, Gate& gate) :
  gate{gate}, syncc{true}, an_size{0}, an_ops{nullptr}, an_bits{nullptr}
{
  qbits.load(path, arma_binary);
  size = log2(qbits.n_rows);
  state = qbits.n_cols > 1 ? "mix" : "pure";
  ops = new char[size];
  bits = new Bit[size]();
  std::memset(ops, 'I', size*sizeof(char));
  std::srand(seed);
}

QSystem::~QSystem() {
  delete[] ops;
  delete[] bits;
  if (an_ops) delete[] an_ops;
  if (an_bits) delete[] an_bits;
}

void QSystem::print_state() {
  auto to_bits = [&](size_t i) {
    std::string sbits{'|'};
    for (size_t j = 0; j < size; j++)
      sbits += i & 1ul << (size+an_size-j-1)? '1' : '0';
    sbits += an_size == 0? ">" : ">|";
    for (size_t j = size; j < size+an_size; j++)
      sbits += i & 1ul << (size+an_size-j-1)? '1' : '0';
    sbits += an_size == 0? "" : ">";
    return sbits;    
  };

  auto cx_to_str = [&](std::complex<double> i) {
    std::stringstream ss;
    if (fabs(i.imag()) < 1e-14) {
      ss << std::showpos << std::fixed
         << std::setprecision(3)  << i.real() << "       ";
    } else if (fabs(i.real()) < 1e-14) {
      ss << std::showpos << std::fixed
         << std::setprecision(3) << std::setw(12) << i.imag() << 'i';
    } else {
      ss << std::showpos << std::fixed << std::setprecision(3) << i.real()
         << i.imag() << 'i';
    }
    return ss.str();
  };

  if (not syncc) sync();
  if (state == "pure") {
    std::stringstream ss;

    for (auto i = qbits.begin(); i != qbits.end(); ++i) {
      if (abs((cx_double)*i) < 1e-14) continue; 
      ss << cx_to_str(*i) << to_bits(i.row()) << '\n';
    }
    std::cout << ss.str() << std::endl; 
  } else if (state == "mix") {
    for (auto i = qbits.begin(); i != qbits.end(); ++i) {
      auto aux = cx_to_str(*i);
      std::cout << "(" << i.row() << ", " << i.col() << ")    " <<
        (aux == ""? "1" : aux)  << std::endl;
    }
  }
}

size_t QSystem::get_size() {
  return size;
}

std::vector<int> QSystem::get_bits() {
  std::vector<int> vec;
  for (size_t i = 0; i < size; i++)
    vec.push_back(bits[i]);
  return vec;
}

size_t QSystem::get_an_size() {
  return an_size;
}

std::vector<int> QSystem::get_an_bits() {
  std::vector<int> vec;
  for (size_t i = 0; i < an_size; i++)
    vec.push_back(an_bits[i]);
  return vec;
}

void QSystem::change_to(std::string state) {
  if (state != "mix" and state != "pure") 
    throw std::invalid_argument{"Argument \'state\' must be \"pure\" or \"mix\", not \""
      + state + "\"."};

  if (state == this->state) 
    return;

  if (state == "mix") {
    qbits = qbits*qbits.t();
  } else if (state == "pure") {
    sp_cx_mat nqbits{1ul << (size+an_size), 1};
    for (size_t i = 0; i < 1ul << (size+an_size); i++)
      nqbits(i,0) = sqrt(qbits(i,i).real());
    qbits = nqbits;
  }
  
  this->state = state;
}

std::string QSystem::get_state() {
  return state;
}

sp_cx_mat QSystem::make_gate(sp_cx_mat gate, size_t qbit) {
  sp_cx_mat m;
  if (qbit == 0) {
    size_t eyesize = 1ul << (size+an_size-1);
    m = kron(gate, eye<sp_cx_mat>(eyesize, eyesize));
  } else if (qbit == size+an_size-1) {
    size_t eyesize = 1ul << (size+an_size-1);
    m = kron(eye<sp_cx_mat>(eyesize, eyesize), gate);
  } else {
    size_t eyesize = 1ul << qbit;
    m = kron(eye<sp_cx_mat>(eyesize, eyesize), gate);
    eyesize = 1ul << (size+an_size-qbit-1);
    m = kron(m, eye<sp_cx_mat>(eyesize, eyesize));
  }
  return m;
}

void QSystem::save(std::string path) {
  if (not syncc) sync();
  qbits.save(path, arma_binary);
}

