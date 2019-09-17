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

#include "../header/using.h"

str utility::cx_to_str(complex i, bool use_sqrt) {
  std::stringstream ss;
  if (fabs(i.imag()) < 1e-14) {
    double tmp = 1.0/pow(i.real(), 2);
    int itmp = tmp;
    if (use_sqrt and itmp != 1 and abs(tmp - double(itmp)) <= 1e-9) {
      ss << (signbit(i.real())? '-' : '+') << std::left << std::setw(11)
         << str{"1/sqrt(" + std::to_string(itmp) + ")"} << std::left
         << std::setw(12) << ' ';
    } else {
      ss << std::showpos << std::fixed
         << std::setprecision(9)  
         << i.real() << std::setw(12) << ' ';
    }
  } else if (fabs(i.real()) < 1e-14) {
    double tmp = 1.0/pow(i.imag(), 2);
    int itmp = tmp;
    if (use_sqrt and itmp != 1 and abs(tmp - double(itmp)) <= 1e-9) {
      ss << std::left << std::setw(12) << ' ' 
         << (signbit(i.imag())? '-' : '+')
         << std::left << std::setw(11) 
         << str{"i/sqrt(" + std::to_string(itmp) + ")"};
    } else{
      ss << std::left << std::setw(12) << ' '  
         << std::showpos << std::fixed
         << std::setprecision(8) << i.imag() << 'i';
    }
  } else {
    double tmp = 1.0/pow(i.real(), 2);
    int itmp = tmp;
    if (use_sqrt and itmp != 1 and abs(tmp - double(itmp)) <= 1e-9) 
      ss << (signbit(i.real())? '-' : '+') 
         << std::left << std::setw(11)
         << str{"1/sqtr(" + std::to_string(itmp) + ")"};
    else
      ss << std::showpos << std::fixed << std::setprecision(9) << i.real();
    tmp = 1.0/pow(i.imag(), 2);
    itmp = tmp;
    if (use_sqrt and itmp != 1 and abs(tmp - double(itmp)) <= 1e-9) 
      ss << (signbit(i.imag())? '-' : '+') 
         << std::left << std::setw(11) 
         << str{"i/sqrt(" + std::to_string(itmp) + ")"};
    else
      ss << std::showpos << std::fixed << std::setprecision(8) 
         << i.imag() << 'i';
  }
  return ss.str();
}

str utility::to_bits(size_t i, size_t qsize, size_t asize) {
  auto size = qsize+asize;
  str sbits{'|'};
  for (size_t j = 0; j < qsize; j++)
    sbits += i & 1ul << (size-j-1)? '1' : '0';
  sbits += asize == 0? ">" : ">|";
  for (size_t j = qsize; j < size; j++)
    sbits += i & 1ul << (size-j-1)? '1' : '0';
  sbits += asize == 0? "" : ">";
  return sbits;    
};

