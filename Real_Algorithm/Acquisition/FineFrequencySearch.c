// FineFrequencySearch.c - Tìm kiếm tần số Doppler chính xác sau khi đã có kết quả Acquisition
#include "FineFrequencySearch.h"
#include "Tools.h"
#include <stdlib.h>
#include <math.h>
#include "Acquisition.h"

AcquisitionResult performFineFrequencySearch(const float* signal_in, const float* local_prn, 
                                             int number_of_samples, float f_sampling, float if_f,
                                             AcquisitionResult coarse_result) 
{
    // 1. Kiểm tra: Nếu bước dò thô không bắt được tín hiệu, trả về nguyên trạng
    if (!coarse_result.is_acquired) {
        return coarse_result;
    }

    // Biến number_of_samples chính là số lượng mẫu trong 1 mili-giây
    int signal_10ms = 10 * number_of_samples; 
    
    // 2. Tính toán offset (độ trễ) để ép khối dữ liệu bắt đầu ở đúng ranh giới mã C/A
    int delay_idx = (number_of_samples - coarse_result.best_code_phase_index) % number_of_samples;

    // Cấp phát bộ nhớ cho mã PRN kéo dài 10ms
    float* long_ca_code = (float*)malloc(signal_10ms * sizeof(float));
    if (!long_ca_code) return coarse_result;

    for (int i = 0; i < signal_10ms; i++) {
        long_ca_code[i] = local_prn[i % number_of_samples];
    }

    // 3. Triệt mã (Wipe-off Code): Nhân dữ liệu thô với mã PRN để chỉ giữ lại sóng mang
    float* x_carrier = (float*)malloc(signal_10ms * sizeof(float));
    if (!x_carrier) { 
        free(long_ca_code); 
        return coarse_result; 
    }

    // Khử nhiễu DC (Zero-mean): Tính giá trị trung bình để triệt tiêu đỉnh nhiễu tại 0 Hz
    float mean = 0.0f;
    for (int i = 0; i < signal_10ms; i++) {
        mean += signal_in[delay_idx + i];
    }
    mean /= signal_10ms;

    // Nhân để triệt mã
    for (int i = 0; i < signal_10ms; i++) {
        x_carrier[i] = (signal_in[delay_idx + i] - mean) * long_ca_code[i];
    }

    // 4. Chuẩn bị FFT: Tính kích thước với kỹ thuật Zero-padding (Gấp 8 lần lũy thừa bậc 2 gần nhất)
    int fft_size = 1;
    while (fft_size < signal_10ms) {
        fft_size *= 2;
    }
    fft_size *= 8;  // Kỹ thuật nội suy phổ mật độ cao, giúp dò Doppler chính xác tới ~1-2 Hz

    Complex* fft_input  = (Complex*)calloc(fft_size, sizeof(Complex));
    Complex* fft_output = (Complex*)malloc(fft_size * sizeof(Complex));
    
    if (!fft_input || !fft_output) {
        free(long_ca_code); 
        free(x_carrier);
        if (fft_input) free(fft_input);
        if (fft_output) free(fft_output);
        return coarse_result;
    }

    // Nạp dữ liệu vào mảng FFT (calloc đã tự động điền 0 vào phần pad)
    for (int i = 0; i < signal_10ms; i++) {
        fft_input[i].real = x_carrier[i];
        fft_input[i].imag = 0.0f; // Tín hiệu thực thì imag = 0
    }

    // Thực thi biến đổi Fourier nhanh (tùy thuộc thư viện bạn đang sử dụng)
    customFFT(fft_input, fft_output, fft_size);

    // =================================================================
    // 5. LỌC PHỔ: Cửa sổ tìm kiếm (Search Window) quanh Coarse Doppler
    // =================================================================
    float target_freq = if_f + coarse_result.best_doppler;
    
    // Tính chỉ số bin trung tâm tương ứng với tần số đích
    int target_bin = (int)roundf(target_freq * fft_size / f_sampling);
    
    // Mở rộng khoảng tìm kiếm +/- 1000 Hz để bao trọn sai số Coarse Search
    int search_radius_bins = (int)ceilf(1000.0f * fft_size / f_sampling); 
    
    int bin_min = target_bin - search_radius_bins;
    int bin_max = target_bin + search_radius_bins;
    
    // Bảo vệ tránh quét tràn mảng FFT (giới hạn từ bin số 1 tới Nyquist)
    if (bin_min < 1) bin_min = 1;
    if (bin_max >= (fft_size / 2)) bin_max = (fft_size / 2) - 1;

    float max_val = 0.0f;
    int max_idx = target_bin; 
    
    // Quét tìm đỉnh năng lượng cao nhất chỉ trong phạm vi cửa sổ
    for (int i = bin_min; i <= bin_max; i++) {
        float mag = fft_output[i].real * fft_output[i].real 
                  + fft_output[i].imag * fft_output[i].imag;
        if (mag > max_val) {
            max_val = mag;
            max_idx = i;
        }
    }

    // 6. Tính toán lại tần số Doppler tinh chỉnh (Fine Doppler)
    float fine_carrier_freq = (float)max_idx * f_sampling / (float)fft_size;
    float fine_doppler = fine_carrier_freq - if_f;

    // Cập nhật kết quả tinh chỉnh vào biến trả về
    coarse_result.best_doppler = fine_doppler;

    // 7. Giải phóng toàn bộ bộ nhớ Heap
    free(long_ca_code); 
    free(x_carrier);
    free(fft_input); 
    free(fft_output);

    // Trả về struct AcquisitionResult đã được cập nhật giá trị Doppler mới
    return coarse_result;
}