// Code Loop Filter - Bộ lọc vòng khóa pha mã (Bậc 2)
#ifndef CODE_LOOP_FILTER_H
#define CODE_LOOP_FILTER_H
#include "CarrierLoopFilter.h" // Cấu trúc LoopFilterState có thể tái sử dụng cho cả PLL và DLL

// Khởi tạo bộ lọc mã 
void CodeLoopFilter_Init(LoopFilterState *state, float b_l, float zeta, float t_int);

// Cập nhật sai số pha mã và trả về lượng hiệu chỉnh tần số mã dư (Hz)
float CodeLoopFilter_Update(LoopFilterState *state, float current_error);

#endif