// Tính toán vị trí bộ thu GPS dựa trên khoảng giả và vị trí vệ tinh
#include "ReceiverPosition.h"
#include <stdio.h>
#include "Tools.h"

#define CLIGHT 299792458.0      // Tốc độ ánh sáng (m/s)
#define F_SAMPLING 5000000.0f   // Tần số lấy mẫu 
#define OMEGA_E_DOT 7.2921151467e-5 // Tốc độ quay của Trái Đất [rad/s]
#include <math.h> // Đảm bảo có math.h để dùng sin, cos
// --- CÁC HÀM BÙ TRỪ KHÍ QUYỂN (ATMOSPHERIC CORRECTIONS) ---

// 1. Mô hình Tầng đối lưu (Troposphere) - Saastamoinen rút gọn
double calculate_tropospheric_delay(double el_rad, double rx_alt) {
    if (el_rad < 0.087) el_rad = 0.087; // Cắt góc ngẩng tối thiểu ở 5 độ để tránh vô cực
    
    // Ở mực nước biển, tín hiệu bị trễ khoảng 2.47 mét tại thiên đỉnh.
    // Càng lên cao, không khí càng loãng -> giảm theo hàm mũ.
    double zenith_delay = 2.47 * exp(-rx_alt / 7000.0); 
    
    // Ánh xạ (Mapping function) theo góc ngẩng của vệ tinh
    return zenith_delay / sin(el_rad);
}

// 2. Mô hình Tầng điện ly (Ionosphere) - Kinh nghiệm (Ban ngày trung bình)
double calculate_ionospheric_delay(double el_rad) {
    if (el_rad < 0.087) el_rad = 0.087;
    // Do chưa giải mã Subframe 4 (Klobuchar), ta dùng mức trễ thiên đỉnh trung bình ~ 5.0m
    double zenith_delay = 5.0; 
    
    // Ánh xạ tầng điện ly (Lớp vỏ ở độ cao 350km)
    return zenith_delay / sin(el_rad); 
}
/**
 * @brief Khởi tạo trạng thái bộ thu GPS
 * @param rx Con trỏ tới cấu trúc ReceiverState cần khởi tạo
 */
void init_receiver(ReceiverState* rx) {
    rx->rx_x = 0.0;
    rx->rx_y = 0.0;
    rx->rx_z = 0.0;
    rx->rx_clock_bias = 0.0;
    rx->num_active_sats = 0;
    for (int i = 0; i < MAX_CHANNELS; i++) {
        rx->channels[i].is_tracked = false;
        rx->channels[i].has_ephemeris = false;
    }
}
/**
 * @brief Tính toán Khoảng giả (Pseudorange) cho tất cả các kênh dựa trên mẫu Preamble
 * @param rx Con trỏ tới trạng thái bộ thu
 * @param preamble_samples Mảng chứa vị trí Sample Index bắt được Preamble của từng kênh
 */
void calculate_pseudoranges(ReceiverState* rx, const uint64_t* preamble_samples) {
    uint64_t min_sample = 0xFFFFFFFFFFFFFFFFULL;
    int closest_sat_idx = -1;

    for (int i = 0; i < MAX_CHANNELS; i++) {
        if (rx->channels[i].is_tracked && rx->channels[i].has_ephemeris) {
            if (preamble_samples[i] < min_sample) {
                min_sample = preamble_samples[i];
                closest_sat_idx = i;
            }
        }
    }
    if (closest_sat_idx == -1) return;

    // Tính toán lấy lại phần lẻ (sub-millisecond) khớp 100% với bản Python
    double min_travel_time_ms = (double)min_sample / (F_SAMPLING / 1000.0);
    double minimum_floor = floor(min_travel_time_ms);

    printf("\n--- TINH TOAN KHOANG GIA (PSEUDORANGE) ---\n");
    printf("Ve tinh neo baseline: PRN %d (Sample: %llu)\n", rx->channels[closest_sat_idx].prn, min_sample);

    for (int i = 0; i < MAX_CHANNELS; i++) {
        if (rx->channels[i].is_tracked && rx->channels[i].has_ephemeris) {
            // Khôi phục công thức Python gốc
            double travel_time_ms = ((double)preamble_samples[i] / (F_SAMPLING / 1000.0)) - minimum_floor + 68.802;
            
            double travel_time_sec = travel_time_ms / 1000.0;
            rx->channels[i].pseudorange = travel_time_sec * CLIGHT;
            
            printf("PRN %2d -> Thoi gian bay: %.6f s | Pseudorange: %11.2f m\n",
                   rx->channels[i].prn, travel_time_sec, rx->channels[i].pseudorange);
        }
    }
    printf("-----------------------------------------\n");
}
/**
 * @brief Tính toán vị trí bộ thu GPS (PVT Solution) dựa trên khoảng giả và vị trí vệ tinh
 * @param rx Con trỏ tới trạng thái bộ thu, chứa thông tin về vệ tinh
 */
int calculate_pvt_solution(ReceiverState* rx) {
    int valid_sats[MAX_CHANNELS];
    int n_sats = 0;
    for (int i = 0; i < MAX_CHANNELS; i++) {
        if (rx->channels[i].is_tracked && rx->channels[i].has_ephemeris) {
            valid_sats[n_sats++] = i;
        }
    }
    if (n_sats < 4) return 0;

    double A[MAX_CHANNELS][4];
    double b[MAX_CHANNELS];
    
    for (int iter = 0; iter < 10; iter++) {
        // -------------------------------------------------------------------
        // Lấy Vĩ độ, Kinh độ, Độ cao hiện tại của bộ thu ở mỗi vòng lặp
        // Để làm tham số đầu vào cho mô hình Khí quyển
        // -------------------------------------------------------------------
        double rx_lat, rx_lon, rx_alt;
        ecef_to_lla(rx->rx_x, rx->rx_y, rx->rx_z, &rx_lat, &rx_lon, &rx_alt);

        for (int i = 0; i < n_sats; i++) {
            SatelliteChannel* sat = &rx->channels[valid_sats[i]];

            // Tịnh tiến Trái Đất (Sagnac effect) - Code cũ của bạn
            double travel_time = sat->pseudorange / CLIGHT;
            double theta = OMEGA_E_DOT * travel_time;
            double sat_x_rot = sat->sat_x * cos(theta) + sat->sat_y * sin(theta);
            double sat_y_rot = -sat->sat_x * sin(theta) + sat->sat_y * cos(theta);
            double sat_z_rot = sat->sat_z;

            double dx = sat_x_rot - rx->rx_x;
            double dy = sat_y_rot - rx->rx_y;
            double dz = sat_z_rot - rx->rx_z;
            double rho_0 = sqrt(dx*dx + dy*dy + dz*dz);

            // ===============================================================
            // TÍNH TOÁN BÙ TRỪ KHÍ QUYỂN (ĐÃ FIX OVERFLOW)
            // ===============================================================
            double tropo_delay = 0.0;
            double iono_delay  = 0.0;
            
            // Ở 2 vòng lặp đầu, tọa độ đang ở Tâm Trái Đất (alt = -6400km)
            // Chỉ kích hoạt mô hình khí quyển từ vòng lặp thứ 3 (iter >= 2) 
            // khi tọa độ đã ngoi lên sát mặt đất để tránh exp() bị vô cực.
            if (iter >= 2) {
                double slon = sin(rx_lon), clon = cos(rx_lon);
                double slat = sin(rx_lat), clat = cos(rx_lat);
                
                // Xoay sang hệ ENU (East, North, Up)
                double E = -slon * dx + clon * dy;
                double N = -slat * clon * dx - slat * slon * dy + clat * dz;
                double U = clat * clon * dx + clat * slon * dy + slat * dz;
                
                double el_rad = atan2(U, sqrt(E*E + N*N));

                // Ép rx_alt luôn >= 0 trước khi đưa vào hàm mũ để an toàn tuyệt đối
                double safe_alt = (rx_alt < 0.0) ? 0.0 : rx_alt; 
                
                tropo_delay = calculate_tropospheric_delay(el_rad, safe_alt);
                iono_delay  = calculate_ionospheric_delay(el_rad);
            }
            // ===============================================================

            A[i][0] = -dx / rho_0;
            A[i][1] = -dy / rho_0;
            A[i][2] = -dz / rho_0;
            A[i][3] = 1.0;

            // Cộng thêm trễ khí quyển (sẽ bằng 0 ở iter 0 và 1)
            double computed_range = rho_0 + rx->rx_clock_bias - (sat->sat_clock_bias * CLIGHT) 
                                  + tropo_delay + iono_delay;
                                  
            b[i] = sat->pseudorange - computed_range;
        }

        double AtA[4][4] = {0};
        double Atb[4] = {0};
        for (int r = 0; r < 4; r++) {
            for (int c = 0; c < 4; c++) {
                for (int k = 0; k < n_sats; k++) AtA[r][c] += A[k][r] * A[k][c];
            }
        }
        for (int r = 0; r < 4; r++) {
            for (int k = 0; k < n_sats; k++) Atb[r] += A[k][r] * b[k];
        }

        double invAtA[4][4];
        if (!invert_matrix_4x4(AtA, invAtA)) return 0;

        double delta_x[4] = {0};
        for (int r = 0; r < 4; r++) {
            for (int c = 0; c < 4; c++) delta_x[r] += invAtA[r][c] * Atb[c];
        }

        rx->rx_x += delta_x[0];
        rx->rx_y += delta_x[1];
        rx->rx_z += delta_x[2];
        rx->rx_clock_bias += delta_x[3];

        if (sqrt(delta_x[0]*delta_x[0] + delta_x[1]*delta_x[1] + delta_x[2]*delta_x[2]) < 1e-3) {
            rx->num_active_sats = n_sats;
            return 1;
        }
    }
    return 0;
}