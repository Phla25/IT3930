// FineFrequencySearch.c
#include "FineFrequencySearch.h"
#include "Tools.h"
#include <stdlib.h>
#include <math.h>
#include "Acquisition.h"

AcquisitionResult performFineFrequencySearch(const Complex* signal_in, const float* local_prn, 
                                             int number_of_samples, float f_sampling, float if_f,
                                             AcquisitionResult coarse_result) 
{
    if (!coarse_result.is_acquired) {
        return coarse_result;
    }

    int signal_10ms = 10 * number_of_samples; 
    int delay_idx = (number_of_samples - coarse_result.best_code_phase_index) % number_of_samples;
    // int delay_idx = coarse_result.best_code_phase_index; // SỬA LẠI ĐỂ TRÁNH NHẦM LẪN VỚI PHASE INDEX TRONG TRACKING

    float* long_ca_code = (float*)malloc(signal_10ms * sizeof(float));
    if (!long_ca_code) return coarse_result;

    for (int i = 0; i < signal_10ms; i++) {
        long_ca_code[i] = local_prn[i % number_of_samples];
    }

    Complex* x_carrier = (Complex*)malloc(signal_10ms * sizeof(Complex));
    if (!x_carrier) { 
        free(long_ca_code); 
        return coarse_result; 
    }

    // Khử nhiễu DC (Zero-mean) ĐỘC LẬP cho cả 2 kênh I và Q
    float mean_i = 0.0f;
    float mean_q = 0.0f;
    for (int i = 0; i < signal_10ms; i++) {
        mean_i += signal_in[delay_idx + i].real;
        mean_q += signal_in[delay_idx + i].imag;
    }
    mean_i /= signal_10ms;
    mean_q /= signal_10ms;

    // Triệt mã PRN trên miền phức
    for (int i = 0; i < signal_10ms; i++) {
        x_carrier[i].real = (signal_in[delay_idx + i].real - mean_i) * long_ca_code[i];
        x_carrier[i].imag = (signal_in[delay_idx + i].imag - mean_q) * long_ca_code[i];
    }

    int fft_size = 1;
    while (fft_size < signal_10ms) {
        fft_size *= 2;
    }
    fft_size *= 8;  

    Complex* fft_input  = (Complex*)calloc(fft_size, sizeof(Complex));
    Complex* fft_output = (Complex*)malloc(fft_size * sizeof(Complex));
    
    if (!fft_input || !fft_output) {
        free(long_ca_code); free(x_carrier);
        if (fft_input) free(fft_input);
        if (fft_output) free(fft_output);
        return coarse_result;
    }

    // Nạp dữ liệu phức vào FFT
    for (int i = 0; i < signal_10ms; i++) {
        fft_input[i].real = x_carrier[i].real;
        fft_input[i].imag = x_carrier[i].imag;
    }

    customFFT(fft_input, fft_output, fft_size);

    // =================================================================
    // LỌC PHỔ: Xử lý cửa sổ tìm kiếm cho phổ tín hiệu phức
    // =================================================================
    float target_freq = if_f + coarse_result.best_doppler;
    
    // Ánh xạ tần số mục tiêu sang Bin. 
    // Công thức: bin = (freq / fs) * N. Xử lý bọc vòng nếu tần số âm.
    int target_bin = (int)roundf(target_freq * fft_size / f_sampling);
    target_bin = (target_bin % fft_size + fft_size) % fft_size; 
    
    int search_radius_bins = (int)ceilf(1000.0f * fft_size / f_sampling); 
    
    float max_val = 0.0f;
    int max_idx = target_bin; 
    
    // Quét tìm đỉnh, áp dụng cơ chế Wrap-around cho mảng FFT
    for (int offset = -search_radius_bins; offset <= search_radius_bins; offset++) {
        // Wrap-around an toàn: bin -1 tương đương với bin (fft_size - 1)
        int current_bin = (target_bin + offset % fft_size + fft_size) % fft_size;
        
        float mag = fft_output[current_bin].real * fft_output[current_bin].real 
                  + fft_output[current_bin].imag * fft_output[current_bin].imag;
                  
        if (mag > max_val) {
            max_val = mag;
            max_idx = current_bin;
        }
    }

    // Giải mã Bin thành Tần số. Nếu Bin > N/2, nó đại diện cho tần số âm.
    float fine_carrier_freq;
    if (max_idx < fft_size / 2) {
        fine_carrier_freq = (float)max_idx * f_sampling / (float)fft_size;
    } else {
        fine_carrier_freq = (float)(max_idx - fft_size) * f_sampling / (float)fft_size;
    }
    
    float fine_doppler = fine_carrier_freq - if_f;
    coarse_result.best_doppler = fine_doppler;

    free(long_ca_code); 
    free(x_carrier);
    free(fft_input); 
    free(fft_output);

    return coarse_result;
}