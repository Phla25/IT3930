// NCO Carrier Generator - Sinh tín hiệu mang dựa trên tần số cục bộ và pha để trộn với tín hiệu đã hạ tần
#include "NCOCarrierGenerator.h"
#include <math.h>

#ifndef PI
#define PI 3.14159265358979323846f
#endif

void NCOCarrierGenerator_Init(NCOCarrierState *state, float init_carrier_freq) {
    state->carrier_freq = init_carrier_freq;
    state->carrier_phase_accummulate = 0.0f; // Khởi đầu pha từ 0
}

void NCOCarrierGenerator_Generate(NCOCarrierState *state, int num_samples, float sampling_freq,
                                  float *cos_carrier, float *sin_carrier) 
{
    float phase_step = 2.0f * PI * state->carrier_freq / sampling_freq;

    for (int i = 0; i < num_samples; i++) {
        // Tính toán sóng mang cục bộ tại mẫu thứ i
        cos_carrier[i] = cosf(state->carrier_phase_accummulate);
        sin_carrier[i] = sinf(state->carrier_phase_accummulate);

        // Cập nhật pha cho mẫu tiếp theo
        state->carrier_phase_accummulate += phase_step;
        if (state->carrier_phase_accummulate >= 2.0f * PI) {
            state->carrier_phase_accummulate -= 2.0f * PI; // Giữ pha trong khoảng [0, 2*PI]
        }
    }
}
