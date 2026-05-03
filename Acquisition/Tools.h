// Các hàm hỗ trợ chung cho các thuật toán Acquisition
#ifndef TOOLS_H
#define TOOLS_H
#include "Acquisition.h"
#define N_FFT 65536 // Kích thước mảng FFT (lũy thừa của 2 gần nhất với 38192)

// Các hàm hỗ trợ
void bitReversal(Complex* data, int n);
void customFFT(const Complex* input, Complex* output, int n);
void customIFFT(const Complex* input, Complex* output, int n);
#endif // TOOLS_H
