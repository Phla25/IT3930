// Carrier Loop Filter - Bộ lọc vòng khóa sóng mang (Bậc 2)
#include "CarrierLoopFilter.h"
#include <math.h>

void CarrierLoopFilter_Init(LoopFilterState *state, float b_l, float zeta, float t_int){
    // Tính tần số tự nhiên (omega_n) dựa trên băng thông nhiễu (b_l) và hệ số cản (zeta)
    // Áp dụng phương trình (7.18)
    float omega_n = (8.0f * zeta * b_l) / (4.0f * zeta * zeta + 1.0f);
    float w_nT = omega_n * t_int;
    
    // Mẫu số chung sinh ra từ phép biến đổi song tuyến tính (Bilinear Transform)
    float denom = 4.0f + 4.0f * zeta * w_nT + w_nT * w_nT; 

    float k_o = t_int; // Độ lợi tích lũy thời gian của NCO sóng mang (0.001)
    float k_d = 1.0f;  // Bộ phân biệt trả về vòng chu kỳ tuyệt đối (1.0)
    
    // Lưu trữ các hệ số c1, c2 vào trạng thái của bộ lọc
    // Áp dụng phương trình (7.16) và (7.17)
    state->c1 = (8.0f * zeta * w_nT) / denom / (k_o * k_d);
    state->c2 = (4.0f * w_nT * w_nT) / denom / (k_o * k_d);
    
    // Khởi tạo trạng thái ban đầu (bộ nhớ trễ z^-1 cho phương trình sai phân)
    state->old_error = 0.0f;
    state->old_output = 0.0f;
}

float CarrierLoopFilter_Update(LoopFilterState *state, float current_error){
    // Thực hiện hàm truyền kỹ thuật số của riêng bộ lọc F(z) = ((c1 + c2) - c1 * z^(-1)) / (1 - z^(-1))
    // Áp dụng phương trình (7.12)
    // Chuyển đổi sang phương trình sai phân miền thời gian:
    // v[n] = v[n-1] + (c1 + c2) * e[n] - c1 * e[n-1]
    
    float current_output = state->old_output 
                         + (state->c1 + state->c2) * current_error 
                         - state->c1 * state->old_error;
                         
    // Cập nhật trạng thái bộ nhớ cho chu kỳ 1ms tiếp theo
    state->old_output = current_output;
    state->old_error = current_error;

    // Trả về lượng hiệu chỉnh tần số Doppler (Hz) để đưa trực tiếp vào Carrier NCO
    return current_output;
}