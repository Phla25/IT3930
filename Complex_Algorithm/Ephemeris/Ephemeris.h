// Tham số vệ tinh được giải mã từ 3 subframe của ephemeris
#ifndef EPHEMERIS_H 
#define EPHEMERIS_H
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>

#define GPS_PI 3.1415926535898
#define MAX_CHANNELS 8   // Số kênh vệ tinh tối đa

// 1. Cấu trúc quản lý thông tin từng kênh vệ tinh đang theo dõi
typedef struct {
    // === THAM SỐ SUBFRAME 1 (Đồng hồ & Trạng thái) ===
    uint32_t weekNumber;
    uint32_t accuracy;
    uint32_t health;
    double   T_GD;
    uint32_t IODC;
    uint32_t t_oc;
    double   a_f2;
    double   a_f1;
    double   a_f0;

    // === THAM SỐ SUBFRAME 2 (Lịch tinh Phần 1) ===
    uint32_t IODE_sf2;
    double   C_rs;
    double   deltan;
    double   M_0;
    double   C_uc;
    double   e;
    double   C_us;
    double   sqrtA;
    uint32_t t_oe;
    
    // === THAM SỐ SUBFRAME 3 (Lịch tinh Phần 2) ===
    double   C_ic;
    double   omega_0;
    double   C_is;
    double   i_0;
    double   C_rc;
    double   omega;
    double   omegaDot;
    uint32_t IODE_sf3;
    double   iDot;

    // === TOW: Time Of Week (giay trong tuan GPS) ===
    // Decode tu HOW word, dung de tinh t_transmit = TOW - pseudorange/c
    double   TOW;

} Ephemeris;

typedef struct {
    int prn;                    // Số hiệu vệ tinh (PRN ID)
    bool is_tracked;            // Đang được khóa tín hiệu hay không
    bool has_ephemeris;         // Đã giải mã đủ 3 Subframe lịch tinh chưa
    
    double pseudorange;         // Khoảng giả (m) - vừa tính ở bước trước
    float  doppler_freq;        // Tần số Doppler (Hz)
    
    Ephemeris eph;              // Bộ lịch tinh giải mã từ Ephemeris.c
    
    // Tọa độ 3D của riêng vệ tinh này
    double sat_x;               
    double sat_y;               
    double sat_z;               
    double sat_clock_bias;      // Sai số đồng hồ vệ tinh
} SatelliteChannel;

void decode_subframe(const uint8_t* subframe, Ephemeris* eph);
#endif // EPHEMERIS_H