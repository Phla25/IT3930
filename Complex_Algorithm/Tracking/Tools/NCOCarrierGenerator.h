// NCO Carrier Generator - Sinh tín hiệu mang dựa trên tần số cục bộ và pha để trộn với tín hiệu đã hạ tần
#ifndef NCO_CARRIER_GENERATOR_H
#define NCO_CARRIER_GENERATOR_H

typedef struct {
    float carrier_freq;      // Tần số sóng mang hiện tại (IF + Doppler)
    float carrier_phase_accummulate; // Bộ tích lũy pha liên tục (từ 0 đến 2*PI)
} NCOCarrierState;

void NCOCarrierGenerator_Init(NCOCarrierState *state, float init_carrier_freq);

// Sinh mảng sóng sin và cos cục bộ để nhân trực tiếp với tín hiệu thô
void NCOCarrierGenerator_Generate(NCOCarrierState *state, int num_samples, float sampling_freq,
                                  float *cos_carrier, float *sin_carrier);

#endif