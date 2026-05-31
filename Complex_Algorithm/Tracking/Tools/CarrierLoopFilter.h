// Carrier Loop Filter - Bộ lọc vòng khóa sóng mang (Bậc 2)
#ifndef CARRIER_LOOP_FILTER_H
#define CARRIER_LOOP_FILTER_H

typedef struct {
    float c1, c2;            // Hệ số bộ lọc số tính từ toán Bilinear/Euler
    float old_error;         // Lưu lại sai số pha của mili-giây trước
    float old_output;       // Lưu lại sai số tần số đã điều chỉnh của mili-giây trước (Doppler lệch)
} LoopFilterState;

// Khởi tạo hệ số bộ lọc dựa trên Băng thông nhiễu (BL), Hệ số cản (Zeta) và Thời gian tích phân (T = 1ms)
void CarrierLoopFilter_Init(LoopFilterState *state, float b_l, float zeta, float t_int);

// Cập nhật sai số mới và trả về giá trị hiệu chỉnh tần số đầu ra (Doppler lệnh)
float CarrierLoopFilter_Update(LoopFilterState *state, float current_error);

#endif