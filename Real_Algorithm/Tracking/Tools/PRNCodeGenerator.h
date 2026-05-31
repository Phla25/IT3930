// PRN Code Generator - Sinh mã C/A dựa trên PRN và tần số lấy mẫu để tạo ra các mảng mã Early, Prompt, Late
#ifndef PRN_CODE_GENERATOR_H
#define PRN_CODE_GENERATOR_H

#include <stdint.h>

typedef struct {
    int prn;
    float code_freq;        // Tần số mã cơ bản (1.023 MHz + Doppler mã)
    float code_phase_accummulate;   // Bộ tích lũy pha mã (Phase Accumulator)
    int32_t code_phase_step; // Bước nhảy pha mã ứng với tần số lấy mẫu
    float ca_code[1023];   // Lưu mảng 1023 chip mã C/A gốc của vệ tinh này
} PRNCodeState;

// Khởi tạo mảng mã C/A dựa trên PRN
void PRNCodeGenerator_Init(PRNCodeState *state, int prn, float sampling_freq);

// Sinh ra 3 mảng mã Early, Prompt, Late tương ứng với số mẫu của 1 block (thường là 1ms)
void PRNCodeGenerator_Generate(PRNCodeState *state, int num_samples, float sampling_freq, 
                               float correlator_spacing, float *early_code, 
                               float *prompt_code, float *late_code);

#endif