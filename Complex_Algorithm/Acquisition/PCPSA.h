// Parallel Code Phase Search Acquisition (PCPSA) algorithm for GPS signal acquisition
#ifndef PCPSA_H
#define PCPSA_H
#include "Acquisition.h"
#include <direct.h> // Đừng quên include dòng này
void save_acq_bin(int prn, int n_fft, Complex* corr_time, int is_first);
AcquisitionResult performPCPSA(
    const Complex* signal_in, // Tín hiệu đầu vào (mảng số thực)
    const float* local_prn, // Mã cục bộ (mảng số thực)
    int number_of_samples, // Số mẫu trong tín hiệu đầu vào
    int number_of_code_phases, // Số code phase cần kiểm tra
    float f_sampling, // Tần số lấy mẫu của tín hiệu đầu vào
    float if_f, // Tần số trung gian của tín hiệu đầu vào
    int prn // PRN code của vệ tinh
);

#endif // PCPSA_H
