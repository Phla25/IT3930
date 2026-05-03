#include <stdio.h>
#include "PCPSA.h"
#include <stdlib.h>
#include <math.h>
#include "Tools.h"
AcquisitionResult performPCPSA(
    const float* signal_in, // Tín hiệu đầu vào (mảng số thực)
    const float* local_prn, // Mã cục bộ (mảng số thực)
    int number_of_samples, // Số mẫu trong tín hiệu đầu vào
    int number_of_code_phases, // Số code phase cần kiểm tra
    float f_sampling, // Tần số lấy mẫu của tín hiệu đầu vào
    float if_f // Tần số trung gian của tín hiệu đầu vào
){
    AcquisitionResult result = {0, 0.0f, 0, 0.0f};
    float max_correlation = 0.0f; 
    float threshold = 500000.0f; 

    Complex* prn_complex = (Complex*)calloc(N_FFT, sizeof(Complex));
    Complex* prn_fft = (Complex*)malloc(N_FFT * sizeof(Complex));
    Complex* signal_baseband = (Complex*)malloc(N_FFT * sizeof(Complex));
    Complex* signal_fft = (Complex*)malloc(N_FFT * sizeof(Complex));
    Complex* cross_corr_freq = (Complex*)malloc(N_FFT * sizeof(Complex));
    Complex* cross_corr_time = (Complex*)malloc(N_FFT * sizeof(Complex));

    if (prn_complex == NULL || prn_fft == NULL || signal_baseband == NULL || signal_fft == NULL || cross_corr_freq == NULL || cross_corr_time == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(prn_complex);
        free(prn_fft);
        free(signal_baseband);
        free(signal_fft);
        free(cross_corr_freq);
        free(cross_corr_time);
        return result;
    }

    // Xử lý mã PRN để tạo thành tín hiệu phức (phần ảo bằng 0)
    for (int i = 0; i < number_of_samples; i++) {
        prn_complex[i].real = local_prn[i];
        prn_complex[i].imag = 0.0f;
    }

    // Zero-padding mã PRN lên kích thước N_FFT
    for (int i = number_of_samples; i < N_FFT; i++) {
        prn_complex[i].real = 0.0f;
        prn_complex[i].imag = 0.0f;
    }

    // Biến đổi FFT của mã PRN
    customFFT(prn_complex, prn_fft, N_FFT);

    // Quét tần số Doppler trong giới hạn +/- 10 kHz
    for (float doppler = -DOPPLER_MAX; doppler <= DOPPLER_MAX; doppler += DOPPLER_STEP) {
        float fc = if_f + doppler;

        for (int i = 0; i < N_FFT; i++) {
            signal_baseband[i].real = 0.0f;
            signal_baseband[i].imag = 0.0f;
        }

        // Hạ tần (Carrier Wipe-off)
        for (int i = 0; i < number_of_samples; i++) {
            float phase = 2.0f * PI * fc * i / f_sampling;
            float i_comp = signal_in[i] * cosf(phase);
            float q_comp = signal_in[i] * (-sinf(phase));
            // Bản copy 1 (nửa đầu)
            signal_baseband[i].real = i_comp;
            signal_baseband[i].imag = q_comp;
            // Bản copy 2 (nửa sau - để hứng đoạn mã PRN bị trượt qua)
            signal_baseband[i + number_of_samples].real = i_comp;
            signal_baseband[i + number_of_samples].imag = q_comp;
        }

        // Nhồi số 0 từ mẫu thứ 38192 đến hết 65535
        for (int i = number_of_samples; i < N_FFT; i++) {
            signal_baseband[i].real = 0.0f;
            signal_baseband[i].imag = 0.0f;
        }

        // Biến đổi FFT của tín hiệu đã hạ tần
        customFFT(signal_baseband, signal_fft, N_FFT);
        // Tính tích chập trong miền tần số (tương đương với nhân phức)
        // Nhân phổ tần số: Tín hiệu * Liên hợp phức của PRN
        for (int i = 0; i < N_FFT; i++) {
            float a = signal_fft[i].real;
            float b = signal_fft[i].imag;
            float c = prn_fft[i].real;
            float d = prn_fft[i].imag; 
            
            cross_corr_freq[i].real = (a * c) + (b * d);
            cross_corr_freq[i].imag = (b * c) - (a * d);
        }

        // Biến đổi ngược IFFT để lấy kết quả tương quan trong miền thời gian
        customIFFT(cross_corr_freq, cross_corr_time, N_FFT);

        // Tìm giá trị tương quan tối đa và vị trí của nó
        for (int i = 0; i < number_of_samples; i++) {
            float i_val = cross_corr_time[i].real;
            float q_val = cross_corr_time[i].imag;
            float corr_power = (i_val * i_val) + (q_val * q_val);
            
            if (corr_power > max_correlation) {
                max_correlation = corr_power;
                result.best_doppler = doppler;
                result.best_code_phase_index = i;
            }
        }
    }
    result.max_correlation = max_correlation;
    if (max_correlation > threshold) {
        result.is_acquired = 1;
    }

    free(prn_complex); free(prn_fft); free(signal_baseband);
    free(signal_fft); free(cross_corr_freq); free(cross_corr_time);

    return result;
}