#include "PRNCodeGenerator.h"
#include <math.h>

static const int G2_SHIFTS[32] = {
    5, 6, 7, 8, 17, 18, 139, 140, 141, 251, 
    252, 254, 255, 256, 257, 258, 469, 470, 471, 472, 
    473, 474, 509, 512, 513, 514, 515, 516, 859, 860, 
    861, 862
};

void PRNCodeGenerator_Init(PRNCodeState *state, int prn, float sampling_freq) {
    state->prn = prn;
    state->code_freq = 1023000.0f;
    state->code_phase_accummulate = 0.0f;

    // FIX 1: Kích thước đúng - thanh ghi 10 bit, index 0-9
    int g1[10], g2[10];
    // FIX 2: Kích thước đúng - lưu 1023 chip thô
    int g1_raw[1023], g2_raw[1023];

    for (int i = 0; i < 10; i++) {
        g1[i] = 1;
        g2[i] = 1;
    }

    for (int i = 0; i < 1023; i++) {
        // FIX 3: Thanh ghi 10 bit → output là bit cuối, index [9]
        g1_raw[i] = g1[9];
        g2_raw[i] = g2[9];

        // FIX 4: Tap đúng theo IS-GPS-200 (đánh số từ 1, nên trừ 1 khi dùng index 0-based)
        // G1: tap 3, 10  → index 2, 9
        int g1_feedback = g1[2] ^ g1[9];
        // G2: tap 2,3,6,8,9,10 → index 1,2,5,7,8,9
        int g2_feedback = g2[1] ^ g2[2] ^ g2[5] ^ g2[7] ^ g2[8] ^ g2[9];

        // Dịch thanh ghi (shift right, đẩy feedback vào đầu)
        for (int j = 9; j > 0; j--) {
            g1[j] = g1[j - 1];
            g2[j] = g2[j - 1];
        }
        g1[0] = g1_feedback;
        g2[0] = g2_feedback;
    }

    // Bước 2: Dịch trễ G2 và XOR với G1
    int shift = G2_SHIFTS[prn - 1];

    for (int i = 0; i < 1023; i++) {
        int g2_delayed_idx = (i - shift + 1023) % 1023;
        int chip = g1_raw[i] ^ g2_raw[g2_delayed_idx];
        state->ca_code[i] = (chip == 0) ? 1.0f : -1.0f;
    }
}

void PRNCodeGenerator_Generate(PRNCodeState *state, int number_of_samples, float f_sampling, 
                               float correlator_spacing, float *early_code, 
                               float *prompt_code, float *late_code) 
{
    // Tính toán bằng float
    float code_step = state->code_freq / (float)f_sampling;
    float spacing = (float)correlator_spacing;

    for (int i = 0; i < number_of_samples; i++) {
        float phase_p = state->code_phase_accummulate;
        float phase_e = phase_p + spacing;
        float phase_l = phase_p - spacing;

        // Xử lý tràn vòng trượt tuần hoàn an toàn
        if (phase_e >= 1023.0f) phase_e -= 1023.0f;
        if (phase_l < 0.0f)     phase_l += 1023.0f;

        // Nội suy zero-order hold bằng cách ép kiểu (int)
        prompt_code[i] = state->ca_code[(int)phase_p];
        early_code[i]  = state->ca_code[(int)phase_e];
        late_code[i]   = state->ca_code[(int)phase_l];

        // Tích lũy pha cho mẫu tiếp theo
        state->code_phase_accummulate += code_step;
        if (state->code_phase_accummulate >= 1023.0f) {
            state->code_phase_accummulate -= 1023.0f;
        }
    }
}