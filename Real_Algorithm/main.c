#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h> 
#include <stdbool.h>
#include <direct.h>
#include <string.h>
#include <pthread.h>

#include "./Acquisition/PCPSA.h"
#include "./Acquisition/FineFrequencySearch.h"
#include "./Tracking/Tracking.h"
#include "./Tracking/Tools/PRNCodeGenerator.h"
#include "./Ephemeris/Preamble.h"
#include "./Ephemeris/Ephemeris.h" 
#include "./Ephemeris/SatellitePosition.h"
#include "./Ephemeris/ReceiverPosition.h"
#include "./Ephemeris/Tools.h"

#define F_SAMPLING 38192000.0f  // fs = 38.192 MHz
#define IF_FREQ     9548000.0f  // if = 9.55 MHz
#define SAMPLES_1MS ((int)(F_SAMPLING * 0.001f))
#define CLIGHT 299792458.0
#define DATA_FILE_PATH "D:\\lessonsatuniversity\\DoAn\\Prj2\\Real_Algorithm\\GPSdata-DiscreteComponents-fs38_192-if9_55.bin"

#define MAX_EPOCHS 10 // Số giây liên tiếp muốn xuất tọa độ liên tục

extern FrameSyncState current_state;
extern int bits_collected;
extern uint8_t word_buffer[62];

pthread_mutex_t preamble_mutex = PTHREAD_MUTEX_INITIALIZER;

// Ma trận toàn cục lưu mốc mẫu tuyệt đối động của từng vệ tinh qua từng giây (Epoch)
uint64_t epoch_samples_matrix[MAX_CHANNELS][MAX_EPOCHS];

typedef struct {
    int prn;
    const float* raw_data;
    int total_samples;
    float init_doppler;
    int initial_phase_idx;
    ReceiverState* rx;
    int ch_idx;
    uint64_t* preamble_sample_idx;
    uint64_t* epoch_row; 
} TrackingThreadArgs;

void track_satellite(int prn, const float* raw_data, int total_samples, float init_doppler, int initial_phase_idx,
                     ReceiverState* rx, int ch_idx, uint64_t* preamble_sample_idx, uint64_t* epoch_row);

void* tracking_worker_routine(void* arg) {
    TrackingThreadArgs* args = (TrackingThreadArgs*)arg;
    track_satellite(
        args->prn, args->raw_data, args->total_samples, 
        args->init_doppler, args->initial_phase_idx, 
        args->rx, args->ch_idx, args->preamble_sample_idx,
        args->epoch_row
    );
    return NULL;
}

int main() {
    _mkdir("Result_Acquisition");
    _mkdir("Result_Tracking");

    printf("==============================================================\n");
    printf(" STARTING RECEIVER GNSS-SDR (PARALLEL MULTI-EPOCH PIPELINE) \n");
    printf("==============================================================\n\n");

    ReceiverState my_receiver;
    init_receiver(&my_receiver);
    uint64_t preamble_samples[MAX_CHANNELS] = {0}; 
    memset(epoch_samples_matrix, 0, sizeof(epoch_samples_matrix));

    FILE* fp = fopen(DATA_FILE_PATH, "rb");
    if (fp == NULL) { return -1; }

    int ms_to_process = 37000; 
    int samples_to_read = SAMPLES_1MS * ms_to_process;
    int8_t* raw_data_8bit = (int8_t*)malloc(samples_to_read * sizeof(int8_t));
    float* raw_signal_buffer = (float*)malloc(samples_to_read * sizeof(float));
    
    fread(raw_data_8bit, sizeof(int8_t), samples_to_read, fp);
    fclose(fp);

    for (int i = 0; i < samples_to_read; i++) raw_signal_buffer[i] = (float)raw_data_8bit[i];
    free(raw_data_8bit); 

    // --- ACQUISITION ---
    int acquired_count = 0;
    int acquired_list[32];
    AcquisitionResult acq_results[32];
    float* local_prn_acq = (float*)malloc(SAMPLES_1MS * sizeof(float));
    float* dummy_buffer = (float*)malloc(SAMPLES_1MS * sizeof(float));

    for (int prn = 1; prn <= 32; prn++) {
        PRNCodeState temp_code_state;
        PRNCodeGenerator_Init(&temp_code_state, prn, F_SAMPLING);
        PRNCodeGenerator_Generate(&temp_code_state, SAMPLES_1MS, F_SAMPLING, 0.0f, dummy_buffer, local_prn_acq, dummy_buffer);
        AcquisitionResult result = performPCPSA(raw_signal_buffer, local_prn_acq, SAMPLES_1MS, SAMPLES_1MS, F_SAMPLING, IF_FREQ, prn);
        result = performFineFrequencySearch(raw_signal_buffer, local_prn_acq, SAMPLES_1MS, F_SAMPLING, IF_FREQ, result);
        if (result.is_acquired) {
            acquired_list[acquired_count] = prn;
            acq_results[acquired_count] = result;
            acquired_count++;
        }
    }
    free(local_prn_acq); free(dummy_buffer);

    if (acquired_count == 0) { free(raw_signal_buffer); return 0; }

    // --- LAUNCH THREAD POOL ---
    int ch_idx = 0;
    pthread_t thread_pool[MAX_CHANNELS];
    TrackingThreadArgs thread_arguments[MAX_CHANNELS];

    for (int i = 0; i < acquired_count; i++) {
        if (ch_idx >= MAX_CHANNELS) break; 
        int prn = acquired_list[i];
        AcquisitionResult res = acq_results[i];

        TrackingThreadArgs* args = &thread_arguments[ch_idx];
        args->prn = prn;
        args->raw_data = raw_signal_buffer;
        args->total_samples = samples_to_read;
        args->init_doppler = res.best_doppler;
        args->initial_phase_idx = res.best_code_phase_index;
        args->rx = &my_receiver;
        args->ch_idx = ch_idx;
        args->preamble_sample_idx = &preamble_samples[ch_idx];
        args->epoch_row = epoch_samples_matrix[ch_idx];

        pthread_create(&thread_pool[ch_idx], NULL, tracking_worker_routine, (void*)args);
        ch_idx++;
    }

    printf("\n>> Processing tracking loops in parallel...\n");
    for (int i = 0; i < ch_idx; i++) {
        pthread_join(thread_pool[i], NULL);
    }

    // ==========================================================
    // BƯỚC 3: GIẢI ĐỊNH VỊ ĐA KỶ NGUYÊN ĐỘNG (ABSOLUTE TIME-TAGGING)
    // ==========================================================
    printf("\n==============================================================\n");
    printf(" START GENERATING CONTINUOUS NAVIGATION POSITION (MULTI-EPOCH) \n");
    printf("==============================================================\n");

    uint64_t temporary_epoch_samples[MAX_CHANNELS];

    for (int epoch = 0; epoch < MAX_EPOCHS; epoch++) {
        for (int ch = 0; ch < MAX_CHANNELS; ch++) {
            temporary_epoch_samples[ch] = epoch_samples_matrix[ch][epoch];
        }

        printf("\n--- TINH TOAN KHOANG GIA (PSEUDORANGE) ---\n");
        calculate_pseudoranges(&my_receiver, temporary_epoch_samples);

        uint64_t min_sample = 0xFFFFFFFFFFFFFFFFULL;
        for (int i = 0; i < MAX_CHANNELS; i++) {
            if (my_receiver.channels[i].is_tracked && my_receiver.channels[i].has_ephemeris) {
                if (temporary_epoch_samples[i] < min_sample) min_sample = temporary_epoch_samples[i];
            }
        }

        // 3. Giải phương trình Kepler dựa trên khoảng dịch mẫu tuyệt đối từ Epoch 0
        for (int i = 0; i < MAX_CHANNELS; i++) {
            if (my_receiver.channels[i].is_tracked && my_receiver.channels[i].has_ephemeris) {
                    
                // 🌟 BẢN VÁ QUYẾT ĐỊNH: Tính thời gian trôi qua thực tế dựa trên hiệu số mẫu của từng kênh
                double dt_seconds = (double)(temporary_epoch_samples[i] - epoch_samples_matrix[i][0]) / F_SAMPLING;
                    
                double time_delay = (double)(temporary_epoch_samples[i] - min_sample) / F_SAMPLING;
                    
                // Thời gian phát chính xác = TOW gốc + Số giây trôi qua thực tế (dt_seconds) - Trễ hình học
                double t_transmit = (double)my_receiver.channels[i].eph.TOW + dt_seconds - time_delay;
                    
                calculate_satellite_position(
                    &my_receiver.channels[i].eph, t_transmit, 
                    &my_receiver.channels[i].sat_x, &my_receiver.channels[i].sat_y, 
                    &my_receiver.channels[i].sat_z, &my_receiver.channels[i].sat_clock_bias
                );
            }
        }

        if (calculate_pvt_solution(&my_receiver)) {
            double final_lat = 0.0, final_lon = 0.0, final_alt = 0.0;
            ecef_to_lla(my_receiver.rx_x, my_receiver.rx_y, my_receiver.rx_z, &final_lat, &final_lon, &final_alt);

            printf("🌟 [EPOCH t0 + %02d sec] -> Vĩ độ: %12.8f° | Kinh độ: %12.8f° | Cao độ: %8.3f m\n", 
                   epoch, final_lat, final_lon, final_alt);
        } else {
            printf("❌ [EPOCH t0 + %02d sec] -> Định vị thất bại.\n", epoch);
        }
    }
    return 0;
}

void track_satellite(int prn, const float* raw_data, int total_samples, float init_doppler, int initial_phase_idx,
                     ReceiverState* rx, int ch_idx, uint64_t* preamble_sample_idx, uint64_t* epoch_row) 
{
    TrackingChannel channel;
    Tracking_Init(&channel, prn, F_SAMPLING, IF_FREQ, init_doppler);

    memset(&rx->channels[ch_idx].eph, 0, sizeof(Ephemeris));
    rx->channels[ch_idx].has_ephemeris = false;

    FrameSyncState local_current_state = SEARCHING_PREAMBLE;
    int local_bits_collected = 0;
    uint8_t local_word_buffer[62];
    memset(local_word_buffer, 0, sizeof(local_word_buffer));

    int delay_idx = (SAMPLES_1MS - initial_phase_idx) % SAMPLES_1MS;
    const float* aligned_signal_stream = raw_data + delay_idx;
    int remaining_samples = total_samples - delay_idx;

    int sample_offset = 0;
    CorrelatorOutputs outputs;
    float current_nav_bit_energy = 0.0f;
    int ms_counter = 0;
    int bit_counter = 0;
    int pull_in_timer = 0;
    uint64_t bit_start_samples[62] = {0};
    uint64_t current_bit_start_sample = 0;

    int  bit_sync_hist[20] = {0};
    bool bit_sync_done  = false;
    int  bit_boundary   = 0;
    float last_1ms_bit  = 0.0f;
    int  sync_ms_count  = 0;
    bool is_aligned     = false;

    bool    preamble_locked    = false;
    bool    invert_subframe    = false;
    uint8_t subframe_bits[300] = {0};
    int     subframe_bit_count = 0;

    bool got_sf1 = false, got_sf2 = false, got_sf3 = false;
    
    // Trả lại bộ đếm ms động thực tế cho mạch DLL bám bắt
    int tracking_ms_counter = 0;
    int current_epoch_idx = 0;

    rx->channels[ch_idx].prn       = prn;
    rx->channels[ch_idx].is_tracked = true;

    while (sample_offset + SAMPLES_1MS <= remaining_samples) {
        uint64_t loop_start_sample = delay_idx + sample_offset;
        float raw_nav_bit_1ms = 0.0f;
        
        Tracking_ProcessBlock(&channel, &aligned_signal_stream[sample_offset],
                              &sample_offset, &outputs, &raw_nav_bit_1ms);

        pull_in_timer++;
        if (pull_in_timer < 100) continue;

        // 🌟 [BẢN VÁ ĐỘNG BÁM BẮT CHÍ MẠNG]:
        // Kích hoạt ngay khi preamble_locked = true. Cứ sau mỗi 1000 chu kỳ bám bắt 1ms, 
        // ta lấy trực tiếp loop_start_sample hiện tại (đã được mạch DLL điều chỉnh theo Code Doppler thực tế).
        if (preamble_locked) {
            tracking_ms_counter++;
            if (tracking_ms_counter % 1000 == 0 && current_epoch_idx < MAX_EPOCHS) {
                epoch_row[current_epoch_idx] = loop_start_sample;
                current_epoch_idx++;
            }
        }

        // ── GIAI ĐOẠN 1: Đồng bộ Bit ──
        if (!bit_sync_done) {
            if (raw_nav_bit_1ms * last_1ms_bit < 0.0f) bit_sync_hist[sync_ms_count % 20]++;
            last_1ms_bit = raw_nav_bit_1ms;
            sync_ms_count++;
            if (sync_ms_count == 1500) {
                int max_val = -1;
                for (int j = 0; j < 20; j++) {
                    if (bit_sync_hist[j] > max_val) { max_val = bit_sync_hist[j]; bit_boundary = j; }
                }
                bit_sync_done = true;
            }
            continue;
        }

        // ── GIAI ĐOẠN 2: Căn chỉnh ranh giới 20ms ──
        if (!is_aligned) {
            if (sync_ms_count % 20 == bit_boundary) {
                is_aligned = true; ms_counter = 1;
                current_nav_bit_energy = raw_nav_bit_1ms;
                current_bit_start_sample = loop_start_sample; 
            }
            sync_ms_count++; continue;
        }

        // ── GIAI ĐOẠN 3: Tích lũy 20ms ──
        current_nav_bit_energy += raw_nav_bit_1ms; ms_counter++; sync_ms_count++;
        if (ms_counter < 20) continue;

        ms_counter = 0; bit_counter++;
        int final_bit_value = (current_nav_bit_energy >= 0.0f) ? 1 : 0;
        current_nav_bit_energy = 0.0f;

        for (int i = 0; i < 61; i++) bit_start_samples[i] = bit_start_samples[i+1];
        bit_start_samples[61] = current_bit_start_sample;
        current_bit_start_sample = delay_idx + sample_offset;

        // ── GIAI ĐOẠN 4: ĐỒNG BỘ KHUNG (MUTEX) ──
        pthread_mutex_lock(&preamble_mutex);
        current_state = local_current_state;
        bits_collected = local_bits_collected;
        memcpy(word_buffer, local_word_buffer, 62);

        process_new_nav_bit(final_bit_value);

        if (current_state == DECODING_SUBFRAME) {
            if (!preamble_locked) {
                preamble_locked = true;
                *preamble_sample_idx = bit_start_samples[2];
                printf("  -> [TIME REFERENCE] Locked Preamble for PRN %d at Sample: %llu\n", prn, *preamble_sample_idx);
                
                // Khóa cứng Epoch 0 tại mốc Preamble vật lý gốc
                epoch_row[0] = bit_start_samples[2];
                current_epoch_idx = 1;
                tracking_ms_counter = 0; // Đặt về 0 để bộ đếm ms chạy đồng bộ từ Preamble gốc đi lên

                uint8_t pre_byte = 0;
                for (int i = 2; i <= 9; i++) pre_byte = (pre_byte << 1) | word_buffer[i];
                invert_subframe = (pre_byte == PREAMBLE_INVERT);
                for (int k = 0; k < 60; k++) {
                    uint8_t b = word_buffer[k + 2];
                    subframe_bits[k] = invert_subframe ? (1 - b) : b;
                }
                subframe_bit_count = 60;
            } 
            else {
                for (int k = 0; k < 61; k++) word_buffer[k] = word_buffer[k + 1];
                word_buffer[61] = final_bit_value;
                uint8_t current_8bit = 0;
                for (int i = 2; i <= 9; i++) current_8bit = (current_8bit << 1) | word_buffer[i];

                if (current_8bit == PREAMBLE_NORMAL || current_8bit == PREAMBLE_INVERT) {
                    if (check_parity(word_buffer) && check_parity(word_buffer + 30)) {
                        bool is_boundary = (subframe_bit_count >= 297 && subframe_bit_count <= 302);
                        if (is_boundary) {
                            invert_subframe = (current_8bit == PREAMBLE_INVERT);
                            while (subframe_bit_count < 300) subframe_bit_count++;
                        }
                    }
                }
                if (subframe_bit_count < 300) {
                    uint8_t b = (uint8_t)final_bit_value;
                    subframe_bits[subframe_bit_count++] = invert_subframe ? (1 - b) : b;
                }
            }
        }
        local_current_state = current_state;
        local_bits_collected = bits_collected;
        memcpy(local_word_buffer, word_buffer, 62);
        pthread_mutex_unlock(&preamble_mutex);

        // Giải mã Subframe
        if (subframe_bit_count == 300 && !rx->channels[ch_idx].has_ephemeris) {
            uint8_t d30_star = 0; 
            for (int w = 0; w < 10; w++) {
                if (d30_star == 1) { for (int i = 0; i < 24; i++) subframe_bits[30 * w + i] ^= 1; }
                d30_star = subframe_bits[30 * w + 29]; 
            }
            int sf_id = (subframe_bits[49] << 2) | (subframe_bits[50] << 1) | subframe_bits[51];
            decode_subframe(subframe_bits, &rx->channels[ch_idx].eph);
            subframe_bit_count = 0;

            if (sf_id == 1 && rx->channels[ch_idx].eph.weekNumber != 0) got_sf1 = true;
            if (sf_id == 2 && rx->channels[ch_idx].eph.sqrtA        != 0) got_sf2 = true;
            if (sf_id == 3)                                                 got_sf3 = true;

            if (got_sf1 && got_sf2 && got_sf3) {
                printf("  [EPHEMERIS READY] PRN %02d locked all subframes. Continuing to track for Epoch positions...\n", prn);
                rx->channels[ch_idx].has_ephemeris = true;
            }
        } else if (subframe_bit_count == 300) {
            subframe_bit_count = 0; 
        }
    }
}