// SatellitePosition.h - Tính toán vị trí vệ tinh dựa trên tham số ephemeris đã giải mã từ 3 subframe
#ifndef SATELLITE_POSITION_H
#define SATELLITE_POSITION_H
#include "Ephemeris.h"
#include <math.h>
void calculate_satellite_position(const Ephemeris* eph, double t, double* satX, double* satY, double* satZ, double* satClockBias);
double calculate_sat_clock_bias(const Ephemeris* eph, double t, double Ek);
#endif // SATELLITE_POSITION_H