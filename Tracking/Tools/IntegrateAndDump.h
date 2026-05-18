// Integrate And Dump - Tích hợp tín hiệu đã hạ tần với các mã PRN để tạo ra các giá trị I và Q, sau đó sử dụng chúng để tính toán độ lệch pha và điều chỉnh tần số cục bộ
#ifndef INTEGRATE_AND_DUMP_H
#define INTEGRATE_AND_DUMP_H

// Định nghĩa cấu trúc lưu 6 giá trị tích phân sau 1ms
typedef struct {
    float I_E, I_P, I_L; // Nhánh đồng pha (Early, Prompt, Late)
    float Q_E, Q_P, Q_L; // Nhánh vuông pha (Early, Prompt, Late)
} CorrelatorOutputs;

// Thực hiện nhân tín hiệu thô với Sóng mang nội và 3 pha mã, sau đó tích phân (cộng dồn)
void IntegrateAndDump_Process(const float *incoming_signal, int num_samples,
                              const float *cos_carrier, const float *sin_carrier,
                              const float *early_code, const float *prompt_code, const float *late_code,
                              CorrelatorOutputs *outputs);

#endif