// Các hàm hỗ trợ chung cho các thuật toán Acquisition
#ifndef TOOLS_H
#define TOOLS_H
#include "Acquisition.h"
#define N_FFT 16384 // Kích thước mảng FFT (lũy thừa của 2 gần nhất với lũy thừa 2 gấn nhất với 5000)

// Các hàm hỗ trợ
void bitReversal(Complex* data, int n);
void customFFT(const Complex* input, Complex* output, int n);
void customIFFT(const Complex* input, Complex* output, int n);
#endif // TOOLS_H
