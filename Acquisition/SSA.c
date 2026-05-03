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
    float threshold = 500000.0f; 

    float* cosine_carrier = (float*)malloc(number_of_samples * sizeof(float)); 
    float* sine_carrier = (float*)malloc(number_of_samples * sizeof(float)); 

    if (cosine_carrier == NULL || sine_carrier == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        if (cosine_carrier) free(cosine_carrier);
        if (sine_carrier) free(sine_carrier);
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
            
            if (correlation_value > max_correlation) { 
                max_correlation = correlation_value; 
                result.best_doppler = doppler; 
                result.best_code_phase_index = phase_shift; 
            }        
        }
    }
    
    result.max_correlation = max_correlation; 
    if (max_correlation > threshold) { 
        result.is_acquired = 1; 
    }
    
    free(cosine_carrier); 
    free(sine_carrier); 
    
    return result;
}