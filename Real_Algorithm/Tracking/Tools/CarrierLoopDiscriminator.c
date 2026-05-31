// Carrier Loop Discriminator - Đo sai số pha giữa sóng mang đầu vào và bản sao sóng mang cục bộ
#include "CarrierLoopDiscriminator.h"
#include <math.h>

#ifndef PI
#define PI 3.14159265358979323846f
#endif

float CarrierLoopDiscriminator_Evaluate(const CorrelatorOutputs *outputs){
    float I_P = outputs->I_P;
    float Q_P = outputs->Q_P;
    // Sử dụng hàm atan2 để tính toán sai số pha sóng mang
    // Kiểm tra chia cho 0 để tránh lỗi
    if (fabs(I_P) < 1e-6){
        return (Q_P >= 0.0f) ? 0.25f : -0.25f; 
        // Nếu I_P gần bằng 0, chỉ dựa vào dấu của Q_P để xác định sai số pha
    }
    float carrier_phase_error = atanf(Q_P/ I_P) / (2.0f * PI);
    return carrier_phase_error;
}