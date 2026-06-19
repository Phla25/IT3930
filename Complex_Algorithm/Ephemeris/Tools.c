// Các hàm hỗ trợ tách triết các tham số vệ tinh từ subframe ephemeris
#include "Tools.h"
// Các hằng số hình học của Trái Đất theo chuẩn WGS-84
#define WGS84_A 6378137.0
#define WGS84_F (1.0 / 298.257223563) // Độ dẹt f
#define PI 3.1415926535898
// Lấy số nguyên KHÔNG DẤU (tương đương bin2dec trong Python)
uint32_t extract_uint(const uint8_t* bits, int start, int length) {
    uint32_t val = 0;
    for(int i = 0; i < length; i++) {
        val = (val << 1) | bits[start + i];
    }
    return val;
}

// Lấy số nguyên CÓ DẤU (tương đương twosComp2dec trong Python)
int32_t extract_int(const uint8_t* bits, int start, int length) {
    uint32_t val = extract_uint(bits, start, length);
    
    // Nếu bit ngoài cùng bên trái là 1 và chiều dài < 32, thực hiện Sign Extension
    if ((bits[start] == 1) && (length < 32)) {
        uint32_t mask = ~0u << length;
        val = val | mask;
    }
    return (int32_t)val;
}

// Hàm bổ trợ: Nghịch đảo ma trận 4x4 bằng phương pháp Gauss-Jordan
int invert_matrix_4x4(double m[4][4], double inv[4][4]) {
    double temp[4][8];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            temp[i][j] = m[i][j];
            temp[i][j + 4] = (i == j) ? 1.0 : 0.0;
        }
    }
    for (int i = 0; i < 4; i++) {
        int max_row = i;
        for (int k = i + 1; k < 4; k++) {
            if (fabs(temp[k][i]) > fabs(temp[max_row][i])) max_row = k;
        }
        if (max_row != i) {
            for (int j = 0; j < 8; j++) {
                double t = temp[i][j];
                temp[i][j] = temp[max_row][j];
                temp[max_row][j] = t;
            }
        }
        if (fabs(temp[i][i]) < 1e-12) return 0;
        double pivot = temp[i][i];
        for (int j = i; j < 8; j++) temp[i][j] /= pivot;
        for (int k = 0; k < 4; k++) {
            if (k != i) {
                double factor = temp[k][i];
                for (int j = i; j < 8; j++) temp[k][j] -= factor * temp[i][j];
            }
        }
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) inv[i][j] = temp[i][j + 4];
    }
    return 1;
}
/**
 * @brief Chuyển đổi ECEF sang Geodetic bằng phương pháp lặp (Theo chuẩn SGK)
 */
void ecef_to_lla(double X, double Y, double Z, double* lat, double* lon, double* alt) {
    // e^2 = (2 - f)*f
    double e2 = (2.0 - WGS84_F) * WGS84_F; 
    
    // 1. Tính Kinh độ lambda (Trực tiếp, không cần lặp)
    // Directly \lambda = arctan(Y/X)
    double lambda = atan2(Y, X); 

    double p = sqrt(X * X + Y * Y);

    // 2. Khởi tạo giá trị ban đầu để vào vòng lặp
    double h = 0.0; // "starting at h = 0"
    
    // Ước lượng phi ban đầu (nếu h = 0)
    double phi = atan2(Z, p * (1.0 - e2)); 
    
    double N_phi;
    double phi_old;

    // 3. Vòng lặp hội tụ cho phi và h
    for (int i = 0; i < 10; i++) {
        phi_old = phi;

        // Tính Bán kính cong đài chuẩn N_phi
        N_phi = WGS84_A / sqrt(1.0 - e2 * sin(phi) * sin(phi));

        // Cập nhật Vĩ độ phi theo công thức (8.37)
        // term = (1 - (e^2 * N_phi) / (N_phi + h))
        double term = 1.0 - (e2 * N_phi) / (N_phi + h);
        phi = atan2(Z, p * term);

        // Cập nhật Độ cao h theo công thức (8.38)
        if (fabs(cos(phi)) < 1e-6) {
            // Tránh chia cho 0 ở hai cực Trái Đất
            h = (Z / sin(phi)) - N_phi * (1.0 - e2); 
        } else {
            // Bình thường
            h = (p / cos(phi)) - N_phi;
        }

        // Kiểm tra hội tụ: Nếu sự thay đổi của góc phi siêu nhỏ thì dừng
        if (fabs(phi - phi_old) < 1e-12) {
            break;
        }
    }

    // 4. Chuyển Radian sang Độ và xuất kết quả
    *lat = phi * 180.0 / PI;
    *lon = lambda * 180.0 / PI;
    *alt = h;
}