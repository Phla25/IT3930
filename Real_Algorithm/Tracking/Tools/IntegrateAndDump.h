// Integrate And Dump - Tích hợp tín hiệu đã hạ tần với các mã PRN để tạo ra các giá trị I và Q, sau đó sử dụng chúng để tính toán độ lệch pha và điều chỉnh tần số cục bộ
#ifndef INTEGRATE_AND_DUMP_H
#define INTEGRATE_AND_DUMP_H

// Thêm vào cấu trúc CorrelatorOutputs để lưu trữ đầy đủ năng lượng 3 nhánh mã
typedef struct {
    float I_E; float Q_E; // Nhánh Sớm (Early)
    float I_P; float Q_P; // Nhánh Đúng pha (Prompt)
    float I_L; float Q_L; // Nhánh Trễ (Late)
    
    // THÊM: Biên độ hình học (Envelope) tính toán bằng sqrt phục vụ vẽ hình 1:1
    float amp_E;
    float amp_P;
    float amp_L;
} CorrelatorOutputs;

// Thực hiện nhân tín hiệu thô với Sóng mang nội và 3 pha mã, sau đó tích phân (cộng dồn)
void IntegrateAndDump_Process(const float *incoming_signal, int num_samples,
                              const float *cos_carrier, const float *sin_carrier,
                              const float *early_code, const float *prompt_code, const float *late_code,
                              CorrelatorOutputs *outputs);

#endif