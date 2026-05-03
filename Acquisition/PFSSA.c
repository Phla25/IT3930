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
    float threshold = 500000.0f; 

    Complex* fft_input = (Complex*)malloc(N_FFT * sizeof(Complex));
    Complex* fft_output = (Complex*)malloc(N_FFT * sizeof(Complex));

    if (fft_input == NULL || fft_output == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        if (fft_input) free(fft_input);
        if (fft_output) free(fft_output);
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

        // Tìm kiếm giá trị tương quan tối đa trong dải Doppler đã định
        for (int bin = (int)min_bin; bin <= (int)max_bin; bin++){
            // Xử lý wrap-around cho chỉ số bin khi nó vượt quá giới hạn của mảng FFT
            int bin_idx = bin;
            if (bin_idx < 0) {
                bin_idx += N_FFT;
            } else if (bin_idx >= N_FFT) {
                bin_idx -= N_FFT;
            }
            float correlation_value = fft_output[bin_idx].real * fft_output[bin_idx].real + fft_output[bin_idx].imag * fft_output[bin_idx].imag; // Tính năng lượng tương quan tại bin này
            if (correlation_value > max_correlation) { 
                max_correlation = correlation_value; 
                result.best_doppler = bin * frequency_resolution - if_f; // Chuyển bin về tần số Doppler
                result.best_code_phase_index = phase_shift;
            }
        }   
    }
    // Kiểm tra nếu giá trị tương quan tối đa vượt qua ngưỡng thì tín hiệu được coi là đã được tìm thấy
    result.max_correlation = max_correlation;
    if (max_correlation > threshold) {
        result.is_acquired = 1; 
    }
    free (fft_input);
    free (fft_output);
    return result;
}


