#include <stdio.h>
#include "PCPSA.h"
#include <stdlib.h>
#include <math.h>
#include "Tools.h"

// BẮT BUỘC PHẢI CÓ 2 DÒNG NÀY ĐỂ NÂNG CẤP BỘ NHỚ LÊN GẤP ĐÔI (Chống trượt Circular Correlation)
#undef N_FFT 
#define N_FFT 131072 // Kích thước mảng FFT (lũy thừa của 2 gần nhất với 76384)

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

    Complex* prn_complex = (Complex*)calloc(N_FFT, sizeof(Complex));
    Complex* prn_fft = (Complex*)malloc(N_FFT * sizeof(Complex));
    Complex* signal_baseband = (Complex*)malloc(N_FFT * sizeof(Complex));
    Complex* signal_fft = (Complex*)malloc(N_FFT * sizeof(Complex));
    Complex* cross_corr_freq = (Complex*)malloc(N_FFT * sizeof(Complex));
    Complex* cross_corr_time = (Complex*)malloc(N_FFT * sizeof(Complex));
    // Mảng 1D lưu năng lượng lớn nhất của từng Pha mã (bất kể Doppler nào) để tìm Đỉnh 2
    float* best_power_per_phase = (float*)calloc(number_of_samples, sizeof(float));

    if (prn_complex == NULL || prn_fft == NULL || signal_baseband == NULL || signal_fft == NULL || cross_corr_freq == NULL || cross_corr_time == NULL || best_power_per_phase == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        if(prn_complex) free(prn_complex);
        if(prn_fft) free(prn_fft);
        if(signal_baseband) free(signal_baseband);
        if(signal_fft) free(signal_fft);
        if(cross_corr_freq) free(cross_corr_freq);
        if(cross_corr_time) free(cross_corr_time);
        if(best_power_per_phase) free(best_power_per_phase);
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

        // Hạ tần (Carrier Wipe-off) và NHÂN ĐÔI TÍN HIỆU
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

        // Nhồi số 0 từ mẫu thứ number_of_samples*2 đến hết N_FFT
        for (int i = number_of_samples*2; i < N_FFT; i++) {
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
        for (int i = 0; i < number_of_samples; i++) { // biến i là độ trễ code phase (Delay)
            float i_val = cross_corr_time[i].real;
            float q_val = cross_corr_time[i].imag;
            float corr_power = (i_val * i_val) + (q_val * q_val); 
            
            // Chuyển i (Độ Trễ - Delay) về Độ Dịch Pha (Phase Shift) để đồng bộ với SSA
            int actual_phase = (number_of_samples - i) % number_of_samples;

            // Cập nhật giá trị lớn nhất cho Pha mã hiện tại vào mảng 1D
            if (corr_power > best_power_per_phase[actual_phase]) {
                best_power_per_phase[actual_phase] = corr_power;
            }
            
            if (corr_power > max_correlation) {
                max_correlation = corr_power;
                result.best_doppler = doppler;
                result.best_code_phase_index = actual_phase;
            }
        }
    }
    
    // ==========================================================
    // BƯỚC RATIO METRIC: TÌM ĐỈNH THỨ 2 VÀ SO SÁNH
    // ==========================================================
    float second_max = 0.0f;
    int exclusion_zone = 38; // Vùng cấm +- 38 mẫu (tương đương 1 chip C/A)

    for(int i = 0; i < number_of_samples; i++) {
        // Tính khoảng cách vòng tròn giữa vị trí i và Đỉnh 1
        int dist = abs(i - result.best_code_phase_index);
        if (dist > number_of_samples / 2) {
            dist = number_of_samples - dist;
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
    
    // Ngưỡng Tỷ số là 2.0 hoặc 2.5
    if (ratio >= 2.4f) {
        result.is_acquired = 1;
    } else {
        result.is_acquired = 0;
    }

    free(prn_complex); 
    free(prn_fft); 
    free(signal_baseband);
    free(signal_fft); 
    free(cross_corr_freq); 
    free(cross_corr_time);
    free(best_power_per_phase);

    return result;
}