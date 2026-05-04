#include <stdio.h>
#include "SSA.h"
#include <stdlib.h>
#include <math.h>

AcquisitionResult performSSA(
    const float* signal_in, // Tín hiệu đầu vào (mảng SỐ THỰC)
    const float* local_prn, // Mã cục bộ (mảng số thực chứa các giá trị +1.0f và -1.0f)
    int number_of_samples,  // Số mẫu trong tín hiệu đầu vào (ví dụ: 38192)
    int number_of_code_phases, // Số pha mã cần kiểm tra
    float f_sampling,       // Tần số lấy mẫu
    float if_f              // Tần số trung gian
) {
    AcquisitionResult result = {0, 0.0f, 0, 0.0f};
    float max_correlation = 0.0f; 
    // float threshold = 10000000.0f; // ĐÃ BỎ: Không dùng ngưỡng cứng nữa, thay bằng Metric Ratio

    float* cosine_carrier = (float*)malloc(number_of_samples * sizeof(float)); 
    float* sine_carrier = (float*)malloc(number_of_samples * sizeof(float)); 
    // Mảng 1D lưu năng lượng lớn nhất của từng Pha mã (bất kể Doppler nào) để tìm Đỉnh 2
    float* best_power_per_phase = (float*)calloc(number_of_code_phases, sizeof(float));

    if (cosine_carrier == NULL || sine_carrier == NULL || best_power_per_phase == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        if (cosine_carrier) free(cosine_carrier);
        if (sine_carrier) free(sine_carrier);
        if (best_power_per_phase) free(best_power_per_phase);
        return result; 
    }

    for (float doppler = -DOPPLER_MAX; doppler <= DOPPLER_MAX; doppler += DOPPLER_STEP) {
        float fc = if_f + doppler; 
        
        // Sinh sóng mang cục bộ
        for (int n = 0; n < number_of_samples; n++) {
            float phase = 2.0f * PI * fc * n / f_sampling; 
            cosine_carrier[n] = cosf(phase); 
            sine_carrier[n] = sinf(phase);
        }

        // Tìm kiếm qua tất cả pha mã
        for (int phase_shift = 0; phase_shift < number_of_code_phases; phase_shift++){
            float i_sum = 0.0f;
            float q_sum = 0.0f;
            
            // Tính tích phân tương quan (tổng từ 0 đến N-1 của x(n)*h(n+k))
            for (int i = 0; i < number_of_samples; i++){
                
                int prn_index = i + phase_shift;
                if (prn_index >= number_of_samples) {
                    prn_index -= number_of_samples;
                }
                
                // Lấy giá trị PRN (+1 hoặc -1) để lột mã (wipe-off)
                float prn_value = local_prn[prn_index]; 
                
                // Hạ tần và triệt tiêu mã 
                float signal_wiped = signal_in[i] * prn_value; 

                i_sum += signal_wiped * cosine_carrier[i]; 
                q_sum += signal_wiped * sine_carrier[i]; 
            }
            
            // Tính năng lượng tương quan
            float correlation_value = (i_sum * i_sum) + (q_sum * q_sum); 
            
            // Cập nhật giá trị năng lượng lớn nhất cho Pha mã hiện tại vào mảng 1D
            if (correlation_value > best_power_per_phase[phase_shift]) {
                best_power_per_phase[phase_shift] = correlation_value;
            }

            if (correlation_value > max_correlation) { 
                max_correlation = correlation_value; 
                result.best_doppler = doppler; 
                result.best_code_phase_index = phase_shift; 
            }        
        }
    }
    
    // ==========================================================
    // BƯỚC RATIO METRIC: TÌM ĐỈNH THỨ 2 VÀ SO SÁNH
    // ==========================================================
    float second_max = 0.0f;
    int exclusion_zone = 38; // Vùng cấm +- 38 mẫu (tương đương 1 chip C/A)

    for(int i = 0; i < number_of_code_phases; i++) {
        // Tính khoảng cách vòng tròn giữa vị trí i và Đỉnh 1
        int dist = abs(i - result.best_code_phase_index);
        if (dist > number_of_code_phases / 2) {
            dist = number_of_code_phases - dist;
        }
        
        // Nếu vị trí i nằm ngoài Vùng cấm, mới được xét làm Đỉnh 2
        if (dist > exclusion_zone) {
            if (best_power_per_phase[i] > second_max) {
                second_max = best_power_per_phase[i];
            }
        }
    }

    // Tính tỷ số (Chống chia cho 0)
    float ratio = 0.0f;
    if (second_max > 0.0f) {
        ratio = max_correlation / second_max;
    }
    
    result.max_correlation = max_correlation; 
    
    // Ngưỡng Tỷ số 2.0 hoặc 2.5
    if (ratio >= 2.4f) { 
        result.is_acquired = 1; 
    } else {
        result.is_acquired = 0;
    }
    
    free(cosine_carrier); 
    free(sine_carrier); 
    free(best_power_per_phase); 
    
    return result;
}