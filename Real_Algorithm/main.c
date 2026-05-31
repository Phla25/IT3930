#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h> // SỬA: Thêm để sử dụng kiểu int8_t chuẩn hệ thống
#include "./Acquisition/PCPSA.h"
#include "./Acquisition/FineFrequencySearch.h"
#include "./Tracking/Tracking.h"
#include "./Tracking/Tools/PRNCodeGenerator.h"
#include <direct.h>

#define F_SAMPLING 38192000.0f  // fs = 38.192 MHz
#define IF_FREQ     9548000.0f  // if = 9.55 MHz
// Số lượng mẫu tương ứng với 1ms dữ liệu thô
#define SAMPLES_1MS ((int)(F_SAMPLING * 0.001f))
#define DATA_FILE_PATH "D:\\lessonsatuniversity\\DoAn\\Prj2\\Real_Algorithm\\GPSdata-DiscreteComponents-fs38_192-if9_55.bin"

void track_satellite(int prn, const float* raw_data, int total_samples, float init_doppler, int initial_phase_idx);

int main() {
    _mkdir("Result_Acquisition");
    _mkdir("Result_Tracking");
    // ==========================================================
    // BƯỚC 1: NẠP FILE DỮ LIỆU THÔ (.BIN) ĐÚNG ĐỊNH DẠNG 8-BIT
    // ==========================================================
    printf("Opening file: %s\n", DATA_FILE_PATH);
    FILE* fp = fopen(DATA_FILE_PATH, "rb");
    if (fp == NULL) {
        fprintf(stderr, "ERROR: Cannot open data file. Please check the path!\n");
        return -1;
    }

    // Đọc trước khoảng 10s dữ liệu để thực hiện bám sát thử nghiệm
    int ms_to_process = 10000;
    int samples_to_read = SAMPLES_1MS * ms_to_process;
    
    // Cấp phát mảng tạm 8-bit để đọc chuẩn xác từng byte dữ liệu từ file
    int8_t* raw_data_8bit = (int8_t*)malloc(samples_to_read * sizeof(int8_t));
    float* raw_signal_buffer = (float*)malloc(samples_to_read * sizeof(float));
    
    if (raw_data_8bit == NULL || raw_signal_buffer == NULL) {
        fprintf(stderr, "LỖI: Cấp phát bộ nhớ thất bại.\n");
        fclose(fp);
        if (raw_data_8bit) free(raw_data_8bit);
        if (raw_signal_buffer) free(raw_signal_buffer);
        return -1;
    }

    // ĐÚNG BẢN CHẤT: Đọc chính xác độ dài 1-byte int8_t từ file
    int samples_read = (int)fread(raw_data_8bit, sizeof(int8_t), samples_to_read, fp);
    fclose(fp);

    // Ép kiểu chuyển đổi dải năng lượng từ int8_t sang float để thuật toán xử lý
    for (int i = 0; i < samples_read; i++) {
        raw_signal_buffer[i] = (float)raw_data_8bit[i];
    }
    
    // Giải phóng ngay mảng tạm 8-bit để tiết kiệm tài nguyên RAM
    free(raw_data_8bit);
    
    printf("Successfully loaded %d samples of real data (~%d ms) into RAM.\n\n", samples_read, ms_to_process);

    // ==========================================================
    // BƯỚC 2: VÒNG LẶP QUÉT TÌM KIẾM TẤT CẢ 32 VỆ TINH GPS
    // ==========================================================
    printf("==============================================================\n");
    printf("   START ACQUISITION (PRN 1 - 32)     \n");
    printf("==============================================================\n");

    int acquired_count = 0;
    int acquired_list[32];
    AcquisitionResult acq_results[32];

    // Tạo các mảng tạm để sinh mã cục bộ 1ms phục vụ pha Acquisition
    float* local_prn_acq = (float*)malloc(SAMPLES_1MS * sizeof(float));
    float* dummy_buffer = (float*)malloc(SAMPLES_1MS * sizeof(float));

    for (int prn = 1; prn <= 32; prn++) {
        PRNCodeState temp_code_state;
        PRNCodeGenerator_Init(&temp_code_state, prn, F_SAMPLING);
        PRNCodeGenerator_Generate(&temp_code_state, SAMPLES_1MS, F_SAMPLING, 0.0f, 
                                   dummy_buffer, local_prn_acq, dummy_buffer);

        // 1. Thực hiện quét thô bằng thuật toán PCPSA (bước nhảy chẵn 500 Hz)
        AcquisitionResult result = performPCPSA(raw_signal_buffer, local_prn_acq, SAMPLES_1MS, SAMPLES_1MS, F_SAMPLING, IF_FREQ, prn);

        // 2. THÊM DÒNG NÀY: Gọi hàm dò mịn bằng Siêu FFT 10ms để tinh chỉnh Doppler về độ phân giải ~4.5 Hz
        result = performFineFrequencySearch(raw_signal_buffer, local_prn_acq, SAMPLES_1MS, F_SAMPLING, IF_FREQ, result);

        // 3. Kiểm tra kết quả (Lúc này kết quả in ra màn hình sẽ là Doppler mịn chuẩn xác)
        if (result.is_acquired) {
            // Log in ra sẽ tự động đổi từ +2000.0 Hz thành con số mịn +1949.2 Hz cực kỳ đẹp mắt
            printf("[FOUND] SV PRN %02d -> Doppler Mịn: %+8.1f Hz | Code Phase Index: %6d\n", 
                   prn, result.best_doppler, (SAMPLES_1MS - result.best_code_phase_index) % SAMPLES_1MS);
            acquired_list[acquired_count] = prn;
            acq_results[acquired_count] = result;
            acquired_count++;
        }
    }

    free(local_prn_acq);
    free(dummy_buffer);

    printf("--------------------------------------------------------------\n");
    printf("RESULT: Acquired %d active SVs.\n", acquired_count);
    printf("==============================================================\n\n");

    // ==========================================================
    // BƯỚC 3: KÍCH HOẠT VÒNG TRACKING CHO CÁC VỆ TINH ĐÃ TÌM THẤY
    // ==========================================================
    if (acquired_count == 0) {
        printf("No SV to track. End program.\n");
        free(raw_signal_buffer);
        return 0;
    }

    for (int i = 0; i < acquired_count; i++) {
        int prn = acquired_list[i];
        AcquisitionResult res = acq_results[i];

        printf(">>> Start Tracking Channel For SV [%02d] <<<\n", prn);
        track_satellite(prn, raw_signal_buffer, samples_read, res.best_doppler, res.best_code_phase_index);
        printf("--------------------------------------------------------------\n");
    }

    free(raw_signal_buffer);
    return 0;
}

void track_satellite(int prn, const float* raw_data, int total_samples, float init_doppler, int initial_phase_idx) {
    TrackingChannel channel;
    Tracking_Init(&channel, prn, F_SAMPLING, IF_FREQ, init_doppler);

    // ==========================================================
    // KHỞI TẠO FILE CSV XUẤT THÔNG SỐ ĐỐI CHIẾU (THÊM CỘT navBit)
    // ==========================================================
    char csv_path[128];
    // Ghi vào thư mục Result_Tracking với đường dẫn dùng dấu / chuẩn
    snprintf(csv_path, sizeof(csv_path), "Result_Tracking/tracking_results_PRN%02d.csv", prn);
    
    FILE* csv_fp = fopen(csv_path, "w");
    if (csv_fp == NULL) {
        fprintf(stderr, "CẢNH BÁO: Không thể tạo file: %s\n", csv_path);
    } else {
        fprintf(csv_fp, "ms,doppler,codeFreq,I_P,Q_P,amp_E,amp_P,amp_L,pllDiscr,dllDiscr,navBit\n");
    }

    // Bù sai lệch pha ban đầu (Phase Shift) -> Căn dòng mẫu bắt đầu tại điểm rìa 1ms
    int delay_idx = (SAMPLES_1MS - initial_phase_idx) % SAMPLES_1MS;
    const float* aligned_signal_stream = raw_data + delay_idx;
    int remaining_samples = total_samples - delay_idx;

    int sample_offset = 0;
    CorrelatorOutputs outputs;
    float current_nav_bit_energy = 0.0f;
    float doppler_sum = 0.0f;
    int doppler_count = 0;
    int ms_counter = 0;
    int bit_counter = 0;
    
    // Trục hoành thời gian tịnh tiến 1ms liên tục phục vụ đồ thị
    int total_ms_elapsed = 0; 

    // Các biến phục vụ thuật toán đồng bộ Bit (Bit Synchronization)
    int bit_sync_found = 0;
    float prev_nav_bit = 0.0f;
    int pull_in_timer = 0; // Bộ đếm thời gian chờ mạch vòng Costas ổn định
    
    // Biến tạm lưu giá trị bit ghi nhận tức thời mỗi ms

    printf("\n >> Bat dau qua trinh bam sat ve tinh PRN %02d...\n", prn);

    while (sample_offset + SAMPLES_1MS <= remaining_samples) {
        float raw_nav_bit_1ms = 0.0f;

        // Xử lý khối tín hiệu 1ms, cập nhật NCO sóng mang và mã nội
        Tracking_ProcessBlock(&channel, &aligned_signal_stream[sample_offset], &sample_offset, 
                              &outputs, &raw_nav_bit_1ms);

        doppler_sum += channel.current_doppler;
        doppler_count++;

        // ----------------------------------------------------------
        // GIAI ĐOẠN 1: Chờ PLL hội tụ & Dò tìm điểm lật pha mã (Bit Sync)
        // ----------------------------------------------------------
        if (!bit_sync_found) {
            pull_in_timer++;
            
            if (pull_in_timer > 100) {
                if ((raw_nav_bit_1ms > 0.0f && prev_nav_bit < 0.0f) || 
                    (raw_nav_bit_1ms < 0.0f && prev_nav_bit > 0.0f)) 
                {
                    bit_sync_found = 1;
                    ms_counter = 1; 
                    current_nav_bit_energy = raw_nav_bit_1ms;
                    printf(" >> [Bit Sync] Ranh gioi bit duoc tim thay tai thoi diem %d ms.\n", pull_in_timer);
                }
            }
            
            // THỰC HIỆN GHI FILE CHO PHA CHƯA ĐỒNG BỘ
            if (csv_fp != NULL) {
                fprintf(csv_fp, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", 
                        total_ms_elapsed++, channel.carrier_state.carrier_freq - IF_FREQ, 
                        channel.code_state.code_freq, outputs.I_P, outputs.Q_P, 
                        outputs.amp_E, outputs.amp_P, outputs.amp_L, 
                        channel.current_carrier_error, channel.current_code_error, outputs.I_P);
            }
            
            prev_nav_bit = raw_nav_bit_1ms;
            continue; 
        }

        // ----------------------------------------------------------
        // GIAI ĐOẠN 2: Gom năng lượng tích phân đủ 20ms để giải mã
        // ----------------------------------------------------------
        current_nav_bit_energy += raw_nav_bit_1ms;
        ms_counter++;


        // THỰC HIỆN GHI FILE CHO PHA ĐÃ ĐỒNG BỘ (CÓ KÈM THEO BIT 0 HOẶC 1)
        if (csv_fp != NULL) {
            fprintf(csv_fp, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", 
                    total_ms_elapsed++, channel.carrier_state.carrier_freq - IF_FREQ, 
                    channel.code_state.code_freq, outputs.I_P, outputs.Q_P, 
                    outputs.amp_E, outputs.amp_P, outputs.amp_L, 
                    channel.current_carrier_error, channel.current_code_error, outputs.I_P);
        }

        if (ms_counter == 20) {
            bit_counter++;
            int final_bit_value = (current_nav_bit_energy >= 0.0f) ? 1 : 0;
            
            printf("  [Bit %03d] Value: %d | Doppler: %+7.1f Hz | DLL: %+5.3f | PLL: %+5.3f\n", 
                bit_counter, final_bit_value, channel.current_doppler,
                channel.current_code_error, channel.current_carrier_error);
                
            current_nav_bit_energy = 0.0f;
            ms_counter = 0;
        }
    }
    
    if (csv_fp != NULL) fclose(csv_fp);
    
    float avg_doppler = doppler_sum / doppler_count;
    printf("\n  >> PRN %02d | Carrier: %.5e Hz | Doppler: %+.0f Hz | Code Phase: %6d samples\n",
       prn, channel.base_if_freq + avg_doppler, avg_doppler, initial_phase_idx);                   
    printf(" >> End of tracking SV %d. Collected %d bits raw data.\n", prn, bit_counter);
}