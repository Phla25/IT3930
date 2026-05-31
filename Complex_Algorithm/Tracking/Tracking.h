#ifndef TRACKING_H
#define TRACKING_H

#include <stdint.h>
#include "./Tools/PRNCodeGenerator.h"
#include "./Tools/NCOCarrierGenerator.h"
#include "./Tools/CarrierLoopFilter.h"
#include "./Tools/IntegrateAndDump.h"

typedef struct {
    int prn;
    float sampling_freq;
    float base_if_freq;      // Tần số trung gian gốc của file dữ liệu
    
    // Khối quản lý trạng thái của các module con
    PRNCodeState code_state;
    NCOCarrierState carrier_state;
    LoopFilterState carrier_filter;
    LoopFilterState code_filter; // Thường dùng khâu tích phân hỗ trợ riêng cho mã
    
    // Các thông số giám sát thời gian thực
    float current_doppler;
    float current_code_error;
    float current_carrier_error;

    float rem_code_phase;   // Số chip dư chưa xử lý từ chu kỳ trước
    float rem_carr_phase;   // Pha sóng mang dư (radian)
} TrackingChannel;

// Khởi tạo toàn bộ kênh bám bắt cho 1 vệ tinh cụ thể
void Tracking_Init(TrackingChannel *channel, int prn, float fs, float if_freq, float init_doppler);

// Xử lý một block tín hiệu thô 1ms, xuất ra giá trị I_P (Dữ liệu định vị) để giải mã bit
void Tracking_ProcessBlock(TrackingChannel *channel, const Complex *signal_block, int* sample_offset, 
                           CorrelatorOutputs *outputs, float *nav_data_bit);

#endif