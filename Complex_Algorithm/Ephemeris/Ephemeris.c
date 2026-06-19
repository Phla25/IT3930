// Tham số vệ tinh được giải mã từ 3 subframe của ephemeris
#include "Ephemeris.h"
#include "Tools.h"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
void decode_subframe(const uint8_t* subframe, Ephemeris* eph) {
    // Đọc Subframe ID (3 bit từ vị trí 49) 
    uint32_t subframeID = extract_uint(subframe, 49, 3);
    // ==========================================================
    // Trích xuất TOW (Từ bit 30, dài 17 bit)
    // TOW trong bản tin chỉ định thời gian của subframe TIẾP THEO. 
    // Vì 1 subframe dài 6 giây, ta trừ 6 để ra thời gian đầu của subframe NÀY.
    // ==========================================================
    uint32_t current_TOW = extract_uint(subframe, 30, 17) * 6 - 6;
    // CHỈ LƯU TOW CỦA SUBFRAME ĐẦU TIÊN BẮT ĐƯỢC (Khớp với mốc preamble_sample_idx)
    if (eph->TOW == 0) {
        eph->TOW = current_TOW;
    }

    // ==========================================================
    // SUBFRAME 1: Thông tin Đồng hồ (Clock Corrections)
    // ==========================================================
    if (subframeID == 1) {
        eph->weekNumber = extract_uint(subframe, 60, 10) + 1024;
        eph->accuracy   = extract_uint(subframe, 72, 4);
        eph->health     = extract_uint(subframe, 76, 6);
        eph->T_GD = (double)extract_int(subframe, 196, 8) * pow(2.0, -31);

        // IODC bị chẻ làm 2 (2 bit MSB và 8 bit LSB)
        uint32_t iodc_msb = extract_uint(subframe, 82, 2);
        uint32_t iodc_lsb = extract_uint(subframe, 210, 8);
        eph->IODC = (iodc_msb << 8) | iodc_lsb;

        eph->t_oc = extract_uint(subframe, 218, 16) * 16;
        eph->a_f2 = (double)extract_int(subframe, 240, 8) * pow(2.0, -55);
        eph->a_f1 = (double)extract_int(subframe, 248, 16) * pow(2.0, -43);
        eph->a_f0 = (double)extract_int(subframe, 270, 22) * pow(2.0, -31);
        
        printf("[OK] Decoded Subframe 1\n");
    } 
    
    // ==========================================================
    // SUBFRAME 2: Lịch tinh Quỹ đạo - Phần 1
    // ==========================================================
    else if (subframeID == 2) {
        eph->IODE_sf2 = extract_uint(subframe, 60, 8);
        eph->C_rs     = (double)extract_int(subframe, 68, 16) * pow(2.0, -5);
        eph->deltan   = (double)extract_int(subframe, 90, 16) * pow(2.0, -43) * GPS_PI;

        // M_0 bị chẻ đôi (Có dấu)
        uint32_t M0_msb = extract_uint(subframe, 106, 8);
        uint32_t M0_lsb = extract_uint(subframe, 120, 24);
        int32_t M0_int = (int32_t)((M0_msb << 24) | M0_lsb);
        eph->M_0 = (double)M0_int * pow(2.0, -31) * GPS_PI;

        eph->C_uc = (double)extract_int(subframe, 150, 16) * pow(2.0, -29);

        // e bị chẻ đôi (Không dấu)
        uint32_t e_msb = extract_uint(subframe, 166, 8);
        uint32_t e_lsb = extract_uint(subframe, 180, 24);
        uint32_t e_int = (e_msb << 24) | e_lsb;
        eph->e = (double)e_int * pow(2.0, -33);

        eph->C_us = (double)extract_int(subframe, 210, 16) * pow(2.0, -29);

        // sqrtA bị chẻ đôi (Không dấu)
        uint32_t sqrtA_msb = extract_uint(subframe, 226, 8);
        uint32_t sqrtA_lsb = extract_uint(subframe, 240, 24);
        uint32_t sqrtA_int = (sqrtA_msb << 24) | sqrtA_lsb;
        eph->sqrtA = (double)sqrtA_int * pow(2.0, -19);

        eph->t_oe = extract_uint(subframe, 270, 16) * 16;
        
        printf("[OK] Decoded Subframe 2\n");
    }

    // ==========================================================
    // SUBFRAME 3: Lịch tinh Quỹ đạo - Phần 2 (BỔ SUNG)
    // ==========================================================
    else if (subframeID == 3) {
        eph->C_ic = (double)extract_int(subframe, 60, 16) * pow(2.0, -29);

        // omega_0 bị chẻ đôi (Có dấu)
        uint32_t omega0_msb = extract_uint(subframe, 76, 8);
        uint32_t omega0_lsb = extract_uint(subframe, 90, 24);
        int32_t omega0_int = (int32_t)((omega0_msb << 24) | omega0_lsb);
        eph->omega_0 = (double)omega0_int * pow(2.0, -31) * GPS_PI;

        eph->C_is = (double)extract_int(subframe, 120, 16) * pow(2.0, -29);

        // i_0 bị chẻ đôi (Có dấu)
        uint32_t i0_msb = extract_uint(subframe, 136, 8);
        uint32_t i0_lsb = extract_uint(subframe, 150, 24);
        int32_t i0_int = (int32_t)((i0_msb << 24) | i0_lsb);
        eph->i_0 = (double)i0_int * pow(2.0, -31) * GPS_PI;

        eph->C_rc = (double)extract_int(subframe, 180, 16) * pow(2.0, -5);

        // omega bị chẻ đôi (Có dấu)
        uint32_t omega_msb = extract_uint(subframe, 196, 8);
        uint32_t omega_lsb = extract_uint(subframe, 210, 24);
        int32_t omega_int = (int32_t)((omega_msb << 24) | omega_lsb);
        eph->omega = (double)omega_int * pow(2.0, -31) * GPS_PI;

        // omegaDot (24 bit liền khối)
        eph->omegaDot = (double)extract_int(subframe, 240, 24) * pow(2.0, -43) * GPS_PI;

        eph->IODE_sf3 = extract_uint(subframe, 270, 8);

        // iDot (14 bit liền khối)
        eph->iDot = (double)extract_int(subframe, 278, 14) * pow(2.0, -43) * GPS_PI;

        printf("[OK] Decoded Subframe 3\n");
    }
}