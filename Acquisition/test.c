#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h> // Để đo thời gian chạy của mỗi thuật toán

#include "Acquisition.h"
#include "SSA.h"
#include "PFSSA.h"
#include "PCPSA.h"

#define FS 38192000.0f     // Tần số lấy mẫu: 38.192 MHz
#define FIF 9550000.0f     // Tần số IF: 9.55 MHz
#define NUM_SAMPLES 38192  // Số mẫu cho 1ms tín hiệu GPS
#define CHIPPING_RATE 1023000.0f // Tốc độ mã chip C/A của GPS (1.023 MHz)

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

// Hàm hỗ trợ in kết quả chung
void print_result(const char* algo_name, AcquisitionResult result, double time_taken) {
    printf("--- [%s] ---\n", algo_name);
    if (result.is_acquired) {
        printf(" => TRANG THAI: Da bat duoc!\n");
        printf(" => Nang luong dinh: %.2f\n", result.max_correlation);
        printf(" => Doppler (Hz):    %.2f\n", result.best_doppler);
        printf(" => Code Phase:      %d\n", result.best_code_phase_index);
    } else {
        printf(" => TRANG THAI: Khong thay! (Dinh cao nhat: %.2f)\n", result.max_correlation);
    }
    printf(" => Thoi gian chay:  %.4f giay\n\n", time_taken);
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

    int target_sv = 21; 
    printf("========================================\n");
    printf("DANG KHAO SAT VE TINH PRN %d...\n", target_sv);
    printf("========================================\n\n");

    int ca_code[1023];
    generate_ca_code(target_sv, ca_code);

    for (int i = 0; i < NUM_SAMPLES; i++) {
        float t = (float)i / FS;
        int chip_index = (int)(t * CHIPPING_RATE) % 1023;
        local_prn_sampled[i] = (ca_code[chip_index] == 0) ? 1.0f : -1.0f;
    }

    // Các biến dùng để đo thời gian
    clock_t start, end;
    double time_taken;
    AcquisitionResult result;

    // ----------------------------------------------------
    // 1. CHẠY THUẬT TOÁN SERIAL SEARCH (SSA)
    // ----------------------------------------------------
    printf("Bat dau chay Serial Search (Rat cham...)...\n");
    start = clock();
    result = performSSA(signal_in, local_prn_sampled, NUM_SAMPLES, NUM_SAMPLES, FS, FIF);
    end = clock();
    time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
    print_result("SERIAL SEARCH (SSA)", result, time_taken);

    // ----------------------------------------------------
    // 2. CHẠY THUẬT TOÁN SONG SONG TẦN SỐ (PFSSA)
    // ----------------------------------------------------
    printf("Bat dau chay Parallel Frequency Space Search...\n");
    start = clock();
    result = performPFSSA(signal_in, local_prn_sampled, NUM_SAMPLES, NUM_SAMPLES, FS, FIF);
    end = clock();
    time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
    print_result("PARALLEL FREQUENCY SEARCH (PFSSA)", result, time_taken);

    // ----------------------------------------------------
    // 3. CHẠY THUẬT TOÁN SONG SONG PHA MÃ (PCPSA)
    // ----------------------------------------------------
    printf("Bat dau chay Parallel Code Phase Search...\n");
    start = clock();
    // Chú ý: Tham số thứ 4 của PCPSA là số lượng code phases. Mình truyền NUM_SAMPLES (38192) để đồng nhất.
    result = performPCPSA(signal_in, local_prn_sampled, NUM_SAMPLES, NUM_SAMPLES, FS, FIF);
    end = clock();
    time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
    print_result("PARALLEL CODE PHASE SEARCH (PCPSA)", result, time_taken);

    printf("========================================\n");
    printf("HOAN THANH KHAO SAT!\n");

    free(signal_in);
    free(local_prn_sampled);
    free(raw_data);

    return 0;
}