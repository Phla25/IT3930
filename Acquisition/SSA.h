// Serial Search Acquisition (SSA) algorithm for GPS signal acquisition
#ifndef SSA_H
#define SSA_H
#include "Acquisition.h"

AcquisitionResult performSSA(
    const float* signal_in, // Tín hiệu đầu vào (mảng số phức)
    const float* local_prn, // Mã cục bộ (mảng số thực)
    int number_of_samples, // Số mẫu trong tín hiệu đầu vào
    int number_of_code_phases, // Số code phase cần kiểm tra
    float f_sampling, // Tần số lấy mẫu của tín hiệu đầu vào
    float if_f // Tần số trung gian của tín hiệu đầu vào
);

#endif // SSA_H
