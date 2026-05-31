// Code Loop Discriminator - Ước lượng độ lệch thời gian (pha mã) giữa mã PRN đầu vào và bản sao mã cục bộ
// Kết quả được dùng để điều khiển Code NCO, giúp căn chỉnh lại tốc độ sinh mã PRN
#include "CodeLoopDiscriminator.h"
#include <math.h>

float CodeLoopDiscriminator_Evaluate(const CorrelatorOutputs *outputs){
    // Sử dụng thuật toán Normalized Early minus Late Amplitude/Envelope
    // Độ lệch pha mã được tính bằng công thức: [sqrt(I_E^2 + Q_E^2) - sqrt(I_L^2 + Q_L^2)]/[sqrt(I_E^2 + Q_E^2) + sqrt(I_L^2 + Q_L^2)]
    // Tính năng lượng nhánh Early
    float power_E = (outputs->I_E * outputs->I_E) + (outputs->Q_E * outputs->Q_E);
    // Tính năng lượng nhánh Late
    float power_L = (outputs->I_L * outputs->I_L) + (outputs->Q_L * outputs->Q_L);
    // Tính toán tử số và mẫu số để chuẩn hóa lỗi (để tránh chia cho 0)
    // Thay power bằng amplitude:
    float amp_E = sqrtf(power_E);
    float amp_L = sqrtf(power_L);
    float normalization = amp_E + amp_L;
    if (normalization < 1e-6f) return 0.0f;
    return (amp_E - amp_L) / normalization; // Trả về giá trị sai số pha mã đã chuẩn hóa
}