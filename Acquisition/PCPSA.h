// Parallel Code Phase Search Acquisition (PCPSA) algorithm for GPS signal acquisition
#ifndef PCPSA_H
#define PCPSA_H
#include "Acquisition.h"
AcquisitionResult performPCPSA(
    const float* signal_in, // Tín hiệu đầu vào (mảng số thực)
    const float* local_prn, // Mã cục bộ (mảng số thực)
    int number_of_samples, // Số mẫu trong tín hiệu đầu vào
    int number_of_code_phases, // Số code phase cần kiểm tra
    float f_sampling, // Tần số lấy mẫu của tín hiệu đầu vào
    float if_f // Tần số trung gian của tín hiệu đầu vào
);

#endif // PCPSA_H
