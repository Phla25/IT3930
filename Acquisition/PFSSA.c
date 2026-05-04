#include <stdio.h>
#include "PFSSA.h"
#include <stdlib.h>
#include <math.h>
#include "Tools.h"

AcquisitionResult performPFSSA(
    const float* signal_in, // Tín hiệu đầu vào (mảng số thực)
    const float* local_prn, // Mã cục bộ (mảng số thực)
    int number_of_samples, // Số mẫu trong tín hiệu đầu vào
    int number_of_code_phases, // Số code phase cần kiểm tra
    float f_sampling, // Tần số lấy mẫu của tín hiệu đầu vào
    float if_f // Tần số trung gian của tín hiệu đầu vào
) {
    AcquisitionResult result = {0, 0.0f, 0, 0.0f};
    float max_correlation = 0.0f; 

    Complex* fft_input = (Complex*)malloc(N_FFT * sizeof(Complex));
    Complex* fft_output = (Complex*)malloc(N_FFT * sizeof(Complex));
    // Mảng 1D lưu năng lượng lớn nhất của từng Pha mã (bất kể Doppler nào) để tìm Đỉnh 2
    float* best_power_per_phase = (float*)calloc(number_of_code_phases, sizeof(float));

    if (fft_input == NULL || fft_output == NULL || best_power_per_phase == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        if (fft_input) free(fft_input);
        if (fft_output) free(fft_output);
        if (best_power_per_phase) free(best_power_per_phase);
        return result; 
    }

    // Độ phân giải tần số
    float frequency_resolution = f_sampling / (float)N_FFT;

    // Tìm kiếm trong giới hạn Doppler +/- 10 kHz xung quanh tần số IF
    float min_bin = (if_f - DOPPLER_MAX) / frequency_resolution;
    float max_bin = (if_f + DOPPLER_MAX) / frequency_resolution;

    // Lặp qua 1023 pha mã
    for (int phase_shift = 0; phase_shift < number_of_code_phases; phase_shift++){
        // Triệt tiêu mã PRN
        for (int i = 0; i < number_of_samples; i++){
            int prn_index = i + phase_shift;
            if (prn_index >= number_of_samples) {
                prn_index -= number_of_samples;
            }
            // Nhân tín hiệu đầu vào với mã cục bộ để triệt tiêu mã
            fft_input[i].real = signal_in[i] * local_prn[prn_index];
            fft_input[i].imag = 0.0f; // Tín hiệu đầu vào là số thực, phần ảo bằng 0
        }

        // Nhồi số 0 từ mẫu thứ 38192 đến hết 65535
        for (int i = number_of_samples; i < N_FFT; i++) {
            fft_input[i].real = 0.0f;
            fft_input[i].imag = 0.0f;
        }

        // Thực hiện FFT để chuyển sang miền tần số (N_FFT điểm là lũy thừa của 2 gần nhất với số mẫu thực tế)
        customFFT(fft_input, fft_output, N_FFT);

        float max_for_this_phase = 0.0f; // Đỉnh cao nhất của pha mã này (để lưu vào mảng 1D)

        // Tìm kiếm giá trị tương quan tối đa trong dải Doppler đã định
        for (int bin = (int)min_bin; bin <= (int)max_bin; bin++){
            // Xử lý wrap-around cho chỉ số bin khi nó vượt quá giới hạn của mảng FFT
            int bin_idx = bin;
            if (bin_idx < 0) {
                bin_idx += N_FFT;
            } else if (bin_idx >= N_FFT) {
                bin_idx -= N_FFT;
            }
            // Tính năng lượng tương quan tại bin này
            float correlation_value = fft_output[bin_idx].real * fft_output[bin_idx].real + fft_output[bin_idx].imag * fft_output[bin_idx].imag; 
            
            // Cập nhật giá trị lớn nhất cho Pha mã hiện tại vào mảng 1D
            if (correlation_value > max_for_this_phase) {
                max_for_this_phase = correlation_value;
            }

            if (correlation_value > max_correlation) { 
                max_correlation = correlation_value; 
                result.best_doppler = bin * frequency_resolution - if_f; // Chuyển bin về tần số Doppler
                result.best_code_phase_index = phase_shift;
            }
        }   
        // Lưu đỉnh cao nhất của Phase này vào mảng
        best_power_per_phase[phase_shift] = max_for_this_phase;
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
    
    free(fft_input);
    free(fft_output);
    free(best_power_per_phase);
    
    return result;
}