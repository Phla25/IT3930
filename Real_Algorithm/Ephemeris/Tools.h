// Các hàm tiện ích để tách triết các tham số vệ tinh
#ifndef TOOLS_H
#define TOOLS_H
#include <stdio.h>
#include <stdint.h>
#include <math.h>
// Lấy số nguyên KHÔNG DẤU (tương đương bin2dec trong Python)
uint32_t extract_uint(const uint8_t* bits, int start, int length);

// Lấy số nguyên CÓ DẤU (tương đương twosComp2dec trong Python)
int32_t extract_int(const uint8_t* bits, int start, int length);

// Hàm bổ trợ: Nghịch đảo ma trận 4x4 bằng phương pháp Gauss-Jordan
int invert_matrix_4x4(double m[4][4], double inv[4][4]);

// Chuyển đổi từ tọa độ ECEF (X, Y, Z) sang LLA (Latitude, Longitude, Altitude)
void ecef_to_lla(double X, double Y, double Z, double* lat, double* lon, double* alt);
#endif // TOOLS_H