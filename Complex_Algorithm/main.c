#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <direct.h>
#include <string.h>

#include "./Acquisition/PCPSA.h"
#include "./Acquisition/FineFrequencySearch.h"
#include "./Tracking/Tracking.h"
#include "./Tracking/Tools/PRNCodeGenerator.h"
#include "./Ephemeris/Preamble.h"
#include "./Ephemeris/Ephemeris.h"
#include "./Ephemeris/SatellitePosition.h"
#include "./Ephemeris/ReceiverPosition.h"
#include "./Ephemeris/Tools.h"

// =======================================================================
// CẤU HÌNH PHẦN CỨNG (Tín hiệu IQ phức, Zero-IF)
// =======================================================================
#define F_SAMPLING   5000000.0f      
#define IF_FREQ      0.0f            
#define SAMPLES_1MS  5000 
#define CLIGHT       299792458.0
#define DATA_FILE_PATH "D:\\lessonsatuniversity\\DoAn\\Prj2\\Complex_Algorithm\\SVTHAYTUNG\\gnss_l1_2ch_sc16_ch0.dat"

#define TOTAL_MS_TO_TRACK 50000 // Chạy tối đa 50 giây 

// Các biến toàn cục từ Preamble.c 
extern FrameSyncState current_state;
extern int            bits_collected;
extern uint8_t        word_buffer[62];

ReceiverState my_receiver;
uint64_t preamble_samples[MAX_CHANNELS] = {0};

// Cấu trúc quản lý luồng dữ liệu song song cho từng Kênh
typedef struct {
    int prn;
    bool is_active;
    FILE* fp;
    
    // [BẢN VÁ]: Khai báo rõ tầng Tracking xử lý tín hiệu DSP
    TrackingChannel tracking_ch; 
    
    // Đồng bộ Bit
    bool bit_sync_done;
    int bit_sync_hist[20];
    int bit_boundary;
    float prev_bit;
    int sync_ms_count;
    bool is_aligned;

    // Tích lũy 20ms
    int ms_counter;
    float nav_bit_energy;

    // Quản lý Ngữ cảnh Toàn cục (Context Backup)
    FrameSyncState frame_state;
    int local_bits_collected;
    uint8_t local_word_buffer[62];
    bool preamble_locked;
    
    // Khung Ephemeris
    uint8_t subframe_bits[300];
    int subframe_bit_count;
    bool got_sf1, got_sf2, got_sf3;

    // Định vị Sample tuyệt đối
    uint64_t global_sample;
    uint64_t bit_start_samples[62];
    uint64_t current_bit_start;
} ChannelTracker;

ChannelTracker trackers[MAX_CHANNELS];

// =======================================================================
// MAIN PIPELINE
// =======================================================================
int main(void) {
    _mkdir("Result_Acquisition");
    _mkdir("Result_Tracking");

    printf("==============================================================\n");
    printf(" GNSS-SDR MULTIPLEXING CORE (1 CHANNEL = 1 SV) \n");
    printf("==============================================================\n\n");

    init_receiver(&my_receiver);
    memset(trackers, 0, sizeof(trackers));

    // ── BƯỚC 1: ACQUISITION ─────────────────────────────────────────────
    int acq_ms = 20;
    int acq_samples = acq_ms * SAMPLES_1MS;

    FILE* fp = fopen(DATA_FILE_PATH, "rb");
    if (!fp) { printf("[ERROR] Cannot open data file!\n"); return -1; }

    int16_t* raw_acq_sc16 = (int16_t*)malloc(acq_samples * 2 * sizeof(int16_t));
    Complex* acq_complex   = (Complex*)malloc(acq_samples * sizeof(Complex));
    
    printf("Loading %d ms Complex IQ for Acquisition...\n", acq_ms);
    fread(raw_acq_sc16, sizeof(int16_t), acq_samples * 2, fp);
    fclose(fp);

    float mean_I = 0.0f, mean_Q = 0.0f;
    for (int i = 0; i < acq_samples; i++) {
        acq_complex[i].real = (float)raw_acq_sc16[2 * i];
        acq_complex[i].imag = (float)raw_acq_sc16[2 * i + 1];
        mean_I += acq_complex[i].real;
        mean_Q += acq_complex[i].imag;
    }
    mean_I /= acq_samples; mean_Q /= acq_samples;
    for (int i = 0; i < acq_samples; i++) {
        acq_complex[i].real -= mean_I;
        acq_complex[i].imag -= mean_Q;
    }
    free(raw_acq_sc16);

    printf("\n>> SCANNING 32 SATELLITES...\n");
    int active_ch = 0;
    PRNCodeState code_state;
    float local_prn[SAMPLES_1MS];

    for (int prn = 1; prn <= 32 && active_ch < MAX_CHANNELS; prn++) {
        PRNCodeGenerator_Init(&code_state, prn, F_SAMPLING);
        for (int i = 0; i < SAMPLES_1MS; i++) {
            local_prn[i] = code_state.ca_code[(int)code_state.code_phase_accummulate];
            code_state.code_phase_accummulate += (1023000.0f / F_SAMPLING);
            if (code_state.code_phase_accummulate >= 1023.0f) code_state.code_phase_accummulate -= 1023.0f;
        }

        AcquisitionResult res = performPCPSA(acq_complex, local_prn, SAMPLES_1MS, SAMPLES_1MS, F_SAMPLING, IF_FREQ, prn);
        if (!res.is_acquired) continue;

        res = performFineFrequencySearch(acq_complex, local_prn, SAMPLES_1MS, F_SAMPLING, IF_FREQ, res);
        printf("[FOUND] PRN %02d -> Doppler: %+8.1f Hz | Code Phase: %5d\n", prn, res.best_doppler, res.best_code_phase_index);

        ChannelTracker* trk = &trackers[active_ch];
        trk->prn = prn;
        trk->is_active = true;
        trk->frame_state = SEARCHING_PREAMBLE;
        trk->global_sample = res.best_code_phase_index;
        
        trk->fp = fopen(DATA_FILE_PATH, "rb");
        fseek(trk->fp, res.best_code_phase_index * 2 * sizeof(int16_t), SEEK_SET);

        // [BẢN VÁ]: Truyền đúng con trỏ cấu trúc TrackingChannel 
        Tracking_Init(&trk->tracking_ch, prn, F_SAMPLING, IF_FREQ, res.best_doppler);
        
        // Khởi tạo tầng Định vị (SatelliteChannel)
        my_receiver.channels[active_ch].prn = prn;
        my_receiver.channels[active_ch].is_tracked = true;
        my_receiver.channels[active_ch].has_ephemeris = false;

        active_ch++;
    }
    free(acq_complex);

    // ── BƯỚC 2: TIME-MULTIPLEXING TRACKING (SONG SONG) ───────────────────
    printf("\n==============================================================\n");
    printf(" START MULTIPLEXING TRACKING (%d CHANNELS SIMULTANEOUSLY)\n", active_ch);
    printf("==============================================================\n");

    int16_t raw_buf[6000 * 2]; 
    Complex comp_buf[6000];
    CorrelatorOutputs outputs;

    bool all_ephemeris_done = false;
    int ms_count = 0;

    while (!all_ephemeris_done && ms_count < TOTAL_MS_TO_TRACK) {
        all_ephemeris_done = true; 
        
        for (int ch = 0; ch < active_ch; ch++) {
            ChannelTracker* trk = &trackers[ch];
            
            // [BẢN VÁ]: Phân tách rõ ràng 2 tầng DSP và Navigation
            TrackingChannel* trk_ch = &trk->tracking_ch;              // Tầng DSP
            SatelliteChannel* rx_ch = &my_receiver.channels[ch];      // Tầng Navigation
            
            if (!trk->is_active || rx_ch->has_ephemeris) continue;
            
            all_ephemeris_done = false; 
            
            float code_phase_step = trk_ch->code_state.code_freq / F_SAMPLING;
            int blksize = (int)ceilf((1023.0f - trk_ch->rem_code_phase) / code_phase_step);
            
            size_t read_items = fread(raw_buf, sizeof(int16_t), blksize * 2, trk->fp);
            if (read_items < (size_t)(blksize * 2)) {
                printf("\n [WARNING] PRN %02d da het file!\n", trk->prn);
                trk->is_active = false;
                continue;
            }
            
            for(int k = 0; k < blksize; k++) {
                comp_buf[k].real = (float)raw_buf[2*k];
                comp_buf[k].imag = (float)raw_buf[2*k + 1];
            }
            
            float raw_nav_bit = 0.0f;
            int dummy_offset = 0;
            Tracking_ProcessBlock(trk_ch, comp_buf, &dummy_offset, &outputs, &raw_nav_bit);
            
            if (!trk->bit_sync_done) {
                if (trk->sync_ms_count > 200 && raw_nav_bit * trk->prev_bit < 0.0f) {
                    trk->bit_sync_hist[trk->sync_ms_count % 20]++;
                }
                trk->prev_bit = raw_nav_bit;
                
                if (trk->sync_ms_count == 1200) {
                    int max_v = 0;
                    for (int j = 0; j < 20; j++) {
                        if (trk->bit_sync_hist[j] > max_v) {
                            max_v = trk->bit_sync_hist[j];
                            trk->bit_boundary = j;
                        }
                    }
                    trk->bit_sync_done = true;
                    printf("  -> [Bit Sync] PRN %02d khoa tai pha %d/20\n", trk->prn, trk->bit_boundary);
                }
                trk->sync_ms_count++;
                trk->global_sample += blksize;
                continue;
            }
            
            if (!trk->is_aligned) {
                if (trk->sync_ms_count % 20 == trk->bit_boundary) {
                    trk->is_aligned = true;
                    trk->ms_counter = 0;
                    trk->nav_bit_energy = 0.0f;
                }
                trk->sync_ms_count++;
                trk->global_sample += blksize;
                continue;
            }
            
            if (trk->ms_counter == 0) {
                trk->current_bit_start = trk->global_sample;
            }
            
            trk->nav_bit_energy += raw_nav_bit;
            trk->ms_counter++;
            trk->sync_ms_count++;
            trk->global_sample += blksize;
            
            if (trk->ms_counter < 20) continue; 
            
            trk->ms_counter = 0;
            int final_bit = (trk->nav_bit_energy >= 0.0f) ? 1 : 0;
            trk->nav_bit_energy = 0.0f;
            
            for(int j = 0; j < 61; j++) trk->bit_start_samples[j] = trk->bit_start_samples[j+1];
            trk->bit_start_samples[61] = trk->current_bit_start;
            
            current_state = trk->frame_state;
            bits_collected = trk->local_bits_collected;
            memcpy(word_buffer, trk->local_word_buffer, 62);
            
            process_new_nav_bit(final_bit); 
            
            trk->frame_state = current_state;
            trk->local_bits_collected = bits_collected;
            memcpy(trk->local_word_buffer, word_buffer, 62);
            
            if (trk->frame_state == DECODING_SUBFRAME && !trk->preamble_locked) {
                trk->preamble_locked = true;
                preamble_samples[ch] = trk->bit_start_samples[2]; 
                printf("  -> [Frame Sync] PRN %02d khoa Preamble tai sample: %llu\n", trk->prn, preamble_samples[ch]);
            }
            
            if (trk->frame_state == DECODING_SUBFRAME) {
                trk->subframe_bits[trk->subframe_bit_count++] = (uint8_t)final_bit;
                
                if (trk->subframe_bit_count == 300) {
                    uint8_t d30_star = 0;
                    for (int w = 0; w < 10; w++) {
                        if (d30_star == 1) {
                            for (int b = 0; b < 24; b++) trk->subframe_bits[30 * w + 2 + b] ^= 1;
                        }
                        d30_star = trk->subframe_bits[30 * w + 29];
                    }

                    int sf_id = (trk->subframe_bits[49] << 2) | (trk->subframe_bits[50] << 1) | trk->subframe_bits[51];

                    // [BẢN VÁ]: Trỏ đúng vào SatelliteChannel để lưu Ephemeris
                    decode_subframe(trk->subframe_bits, &rx_ch->eph);
                    trk->subframe_bit_count = 0;

                    if (sf_id == 1 && rx_ch->eph.weekNumber != 0) trk->got_sf1 = true;
                    if (sf_id == 2 && rx_ch->eph.sqrtA        != 0) trk->got_sf2 = true;
                    if (sf_id == 3)                                 trk->got_sf3 = true;

                    if (trk->got_sf1 && trk->got_sf2 && trk->got_sf3) {
                        printf("\n >>> OK! PRN %02d HOAN THANH EPHEMERIS! <<<\n", trk->prn);
                        rx_ch->has_ephemeris = true; // [BẢN VÁ]
                        fclose(trk->fp); 
                    }
                }
            }
        } 
        
        ms_count++;
        if (ms_count % 5000 == 0 && !all_ephemeris_done) {
            printf("."); fflush(stdout); 
        }
    }

    // ── BƯỚC 3: PVT LEAST SQUARES (TÌM TÒA B1) ─────────────────────────────
    printf("\n\n==============================================================\n");
    printf(" CALCULATING RECEIVER POSITION (PVT)\n");
    printf("==============================================================\n");

    calculate_pseudoranges(&my_receiver, preamble_samples);

    printf("\n--- TINH VI TRI VE TINH ---\n");
    for (int i = 0; i < MAX_CHANNELS; i++) {
        if (!my_receiver.channels[i].is_tracked || !my_receiver.channels[i].has_ephemeris) continue;

        double PR          = my_receiver.channels[i].pseudorange;
        double TOW         = (double)my_receiver.channels[i].eph.TOW;
        double t_transmit  = TOW - PR / CLIGHT;

        calculate_satellite_position(&my_receiver.channels[i].eph, t_transmit,
            &my_receiver.channels[i].sat_x, &my_receiver.channels[i].sat_y,
            &my_receiver.channels[i].sat_z, &my_receiver.channels[i].sat_clock_bias);
            
        printf("PRN %02d: TOW=%.0f | sat_clk=%.6e s\n", 
               my_receiver.channels[i].prn, TOW, my_receiver.channels[i].sat_clock_bias);
    }

    for (int i = 0; i < MAX_CHANNELS; i++) {
        if (!my_receiver.channels[i].is_tracked || !my_receiver.channels[i].has_ephemeris) continue;

        double PR          = my_receiver.channels[i].pseudorange;
        double TOW         = (double)my_receiver.channels[i].eph.TOW;
        double corrected_PR = PR + my_receiver.channels[i].sat_clock_bias * CLIGHT;
        double t_transmit   = TOW - corrected_PR / CLIGHT;

        calculate_satellite_position(&my_receiver.channels[i].eph, t_transmit,
            &my_receiver.channels[i].sat_x, &my_receiver.channels[i].sat_y,
            &my_receiver.channels[i].sat_z, &my_receiver.channels[i].sat_clock_bias);
    }

    if (calculate_pvt_solution(&my_receiver)) {
        double lat, lon, alt;
        ecef_to_lla(my_receiver.rx_x, my_receiver.rx_y, my_receiver.rx_z, &lat, &lon, &alt);

        printf("\n🌟 POSITIONING SUCCESSFUL! 🌟\n");
        printf("  ECEF X: %15.3f m\n", my_receiver.rx_x);
        printf("  ECEF Y: %15.3f m\n", my_receiver.rx_y);
        printf("  ECEF Z: %15.3f m\n", my_receiver.rx_z);
        printf("  Clock Error: %11.9f s\n", my_receiver.rx_clock_bias / CLIGHT);
        printf("\n  🌍 WGS-84 MAP COORDINATES:\n");
        printf("   Latitude  (Vĩ độ) : %14.8f deg\n", lat);
        printf("   Longitude (Kinh độ): %14.8f deg\n", lon);
        printf("   Altitude  (Độ cao) : %14.3f m\n\n", alt);
        
        printf("Bạn có thể copy Vĩ độ và Kinh độ này dán vào Google Maps để xem có đúng Tòa nhà B1 - SoICT không nhé!\n");
    } else {
        printf("\n[FAILED] Cần ít nhất 4 vệ tinh thu thập đủ dữ liệu Ephemeris để định vị!\n");
    }

    return 0;
}