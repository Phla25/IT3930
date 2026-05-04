#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h> 

#include "Acquisition.h"
#include "SSA.h"
#include "PFSSA.h"
#include "PCPSA.h"

#define FS 38192000.0f     
#define FIF 9548000.0f     // Tần số IF CHUẨN: 9.548 MHz
#define NUM_SAMPLES 38192  
#define CHIPPING_RATE 1023000.0f 

const int G2_TAPS[32][2] = {
    {2,6}, {3,7}, {4,8}, {5,9}, {1,9}, {2,10}, {1,8}, {2,9}, {3,10}, {2,3}, 
    {3,4}, {4,5}, {5,6}, {6,7}, {7,8}, {8,9}, {9,10}, {1,4}, {2,5}, {3,6}, 
    {4,7}, {5,8}, {6,9}, {1,3}, {2,4}, {3,5}, {4,6}, {5,7}, {6,8}, {7,9}, 
    {8,10}, {4,10}
};

void generate_ca_code(int sv_id, int* ca_code){
    int g1[10] = {1,1,1,1,1,1,1,1,1,1};
    int g2[10] = {1,1,1,1,1,1,1,1,1,1};

    int t1 = G2_TAPS[sv_id - 1][0] - 1;
    int t2 = G2_TAPS[sv_id - 1][1] - 1;
    for (int i = 0; i < 1023; i++) {
        int g2_out = g2[t1] ^ g2[t2]; 
        ca_code[i] = g1[9] ^ g2_out;
        int fb1 = g1[2] ^ g1[9];
        int fb2 = g2[1] ^ g2[2] ^ g2[5] ^ g2[7] ^ g2[8] ^ g2[9];

        for (int j = 9; j > 0; j--) {
            g1[j] = g1[j-1];
            g2[j] = g2[j-1];
        }
        g1[0] = fb1;
        g2[0] = fb2;
    }
}

// Hàm in kết quả thu gọn để log không bị rối
void print_short_result(const char* algo_name, AcquisitionResult result, double time_taken) {
    if (result.is_acquired) {
        printf("   [%s] THÀNH CÔNG! Đỉnh: %.0f | Doppler: %.0f Hz | Code Phase: %d | Time: %.2fs\n", 
               algo_name, result.max_correlation, result.best_doppler, result.best_code_phase_index, time_taken);
    } else {
        printf("   [%s] Không thấy (Đỉnh nhiễu: %.0f) | Time: %.2fs\n", algo_name, result.max_correlation, time_taken);
    }
}

int main() {
    float* signal_in = (float*)malloc(NUM_SAMPLES * sizeof(float));
    float* local_prn_sampled = (float*)malloc(NUM_SAMPLES * sizeof(float));
    int8_t* raw_data = (int8_t*)malloc(NUM_SAMPLES * sizeof(int8_t)); 

    if (!signal_in || !local_prn_sampled || !raw_data) {
        printf("Lỗi cấp phát bộ nhớ!\n");
        return -1;
    }

    const char* filename = "D:\\lessonsatuniversity\\DoAn\\Prj2\\GPSdata-DiscreteComponents-fs38_192-if9_55.bin"; 
    FILE* file = fopen(filename, "rb");
    if (!file) {
        printf("Khong mo duoc file %s\n", filename);
        return -1;
    }
    
    fread(raw_data, sizeof(int8_t), NUM_SAMPLES, file);
    fclose(file);
    
    for (int i = 0; i < NUM_SAMPLES; i++) {
        signal_in[i] = (float)raw_data[i];
    }

    printf("=================================================================\n");
    printf(" QUET 32 VE TINH BANG CA 3 THUAT TOAN!\n");
    printf("=================================================================\n\n");

    clock_t total_start = clock();
    int found_ssa = 0, found_pfssa = 0, found_pcpsa = 0;

    // VÒNG LẶP DUYỆT 32 VỆ TINH
    for (int sv = 1; sv <= 32; sv++) {
        printf(">>> ĐANG QUÉT VỆ TINH PRN %02d...\n", sv);

        int ca_code[1023];
        generate_ca_code(sv, ca_code);

        for (int i = 0; i < NUM_SAMPLES; i++) {
            float t = (float)i / FS;
            int chip_index = (int)(t * CHIPPING_RATE) % 1023;
            local_prn_sampled[i] = (ca_code[chip_index] == 0) ? 1.0f : -1.0f;
        }

        clock_t start, end;
        double time_taken;
        AcquisitionResult result;

        // 1. SSA
        start = clock();
        result = performSSA(signal_in, local_prn_sampled, NUM_SAMPLES, NUM_SAMPLES, FS, FIF);
        end = clock();
        time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
        print_short_result("SSA  ", result, time_taken);
        if (result.is_acquired) found_ssa++;

        // 2. PFSSA
        // start = clock();
        // result = performPFSSA(signal_in, local_prn_sampled, NUM_SAMPLES, NUM_SAMPLES, FS, FIF);
        //end = clock();
        // time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
        // print_short_result("PFSSA", result, time_taken);
        // if (result.is_acquired) found_pfssa++;

        // 3. PCPSA
        start = clock();
        result = performPCPSA(signal_in, local_prn_sampled, NUM_SAMPLES, NUM_SAMPLES, FS, FIF);
        end = clock();
        time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
        print_short_result("PCPSA", result, time_taken);
        if (result.is_acquired) found_pcpsa++;

        printf("-----------------------------------------------------------------\n");
    }

    clock_t total_end = clock();
    double total_time = ((double) (total_end - total_start)) / CLOCKS_PER_SEC;

    printf("\n=================================================================\n");
    printf(" TONG KET:\n");
    printf(" - SSA tim thay:   %d ve tinh\n", found_ssa);
    // printf(" - PFSSA tim thay: %d ve tinh\n", found_pfssa);
    printf(" - PCPSA tim thay: %d ve tinh\n", found_pcpsa);
    printf("\n TONG THOI GIAN CHAY: %.2f giay (Khang %.2f gio)\n", total_time, total_time / 3600.0);
    printf("=================================================================\n");

    free(signal_in);
    free(local_prn_sampled);
    free(raw_data);

    return 0;
}