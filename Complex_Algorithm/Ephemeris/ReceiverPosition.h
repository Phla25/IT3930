// Tính toán vị trí bộ thu GPS dựa trên khoảng giả và vị trí vệ tinh
#ifndef RECEIVER_POSITION_H
#define RECEIVER_POSITION_H
#include "SatellitePosition.h"

typedef struct {
    // Tọa độ ECEF của chính người dùng (Bộ thu) cần tìm
    double rx_x;                
    double rx_y;                
    double rx_z;                
    
    double rx_clock_bias;       // Sai số đồng hồ bộ thu (quy ra mét hoặc giây)
    double receiver_time;       // Thời gian hiện tại của bộ thu (TOW)
    uint32_t gps_week;          // Số tuần GPS
    
    // Mảng chứa 8 kênh vệ tinh đang chạy song song
    SatelliteChannel channels[MAX_CHANNELS];
    int num_active_sats;        // Số vệ tinh đang hoạt động thực tế
} ReceiverState;

// Khai báo các hàm liên quan để tầng khác gọi
void init_receiver(ReceiverState* rx);
void calculate_pseudoranges(ReceiverState* rx, const uint64_t* preamble_samples);
int calculate_pvt_solution(ReceiverState* rx);

#endif // RECEIVER_H