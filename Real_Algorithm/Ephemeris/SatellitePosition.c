// SatellitePosition.c - Tính toán vị trí vệ tinh dựa trên tham số ephemeris đã giải mã từ 3 subframe
#include "SatellitePosition.h"

#define MU_EARTH 3.986005e14       // Hằng số hấp dẫn của Trái Đất (GM) [m^3/s^2]
#define OMEGA_E_DOT 7.2921151467e-5 // Tốc độ quay của Trái Đất [rad/s]

void calculate_satellite_position(const Ephemeris* eph, double t, double* satX, double* satY, double* satZ, double* satClockBias) {
    // Tính thời gian trôi qua kể từ mốc thời gian ephemeris (t_oe)
    double tk = t - eph->t_oe;
    if (tk > 302400) tk -= 604800; // Nếu tk > 3.5 ngày, trừ đi một tuần
    else if (tk < -302400) tk += 604800; // Nếu tk < -3.5 ngày, cộng thêm một tuần
    
    // Giải phương trình Kepler tìm điểm dị thường
    double a = eph->sqrtA * eph->sqrtA; // Bán trục lớn
    double n0 = sqrt(MU_EARTH / (a * a * a)); //  Vận tốc góc trung bình lý thuyết
    double n = n0 + eph->deltan; // Vận tốc góc trung bình đã hiệu chỉnh
    double Mk = eph->M_0 + n * tk;   // Dị thường trung bình

    // Giải phương trình Kepler bằng phương pháp Newton-Raphson
    double Ek = Mk; // Khởi tạo dị thường lồi bằng dị thường trung bình
    double Ek_prev;
    int iter = 0;
    do {
        Ek_prev = Ek;
        Ek = Mk + eph->e * sin(Ek); // Cập nhật dị thường lồi
        iter++;
    } while (fabs(Ek - Ek_prev) > 1e-12 && iter < 15); // Lặp đến khi hội tụ hoặc đạt max iter

    double sin_fk = (sqrt(1 - eph->e * eph->e) * sin(Ek)) / (1 - eph->e * cos(Ek));
    double cos_fk = (cos(Ek) - eph->e) / (1 - eph->e * cos(Ek));
    double fk = atan2(sin_fk, cos_fk); // Dị thường thật

    double Phi_k = fk + eph->omega; // Góc phương vị
    double uk = Phi_k + eph->C_uc * cos(2 * Phi_k) + eph->C_us * sin(2 * Phi_k); // Góc phương vị đã hiệu chỉnh
    double rk = a * (1 - eph->e * cos(Ek)) + eph->C_rc * cos(2 * Phi_k) + eph->C_rs * sin(2 * Phi_k); // Bán kính đã hiệu chỉnh
    double ik = eph->i_0 + eph->iDot * tk + eph->C_ic * cos(2 * Phi_k) + eph->C_is * sin(2 * Phi_k); // Góc nghiêng đã hiệu chỉnh

    // Vị trí trong mặt phẳng quỹ đạo
    double xk_prime = rk * cos(uk);
    double yk_prime = rk * sin(uk);

    // Kinh độ Điểm nút lên (Longitude of Ascending Node)
    // Trừ đi OMEGA_E_DOT để tính đến việc Trái Đất đang quay bên dưới vệ tinh
    double Omega_k = eph->omega_0 + (eph->omegaDot - OMEGA_E_DOT) * tk - OMEGA_E_DOT * eph->t_oe;

    // Xoay mặt phẳng quỹ đạo vào hệ tọa độ 3D của Trái Đất (ECEF)
    *satX = xk_prime * cos(Omega_k) - yk_prime * cos(ik) * sin(Omega_k);
    *satY = xk_prime * sin(Omega_k) + yk_prime * cos(ik) * cos(Omega_k);
    *satZ = yk_prime * sin(ik);
    *satClockBias = calculate_sat_clock_bias(eph, t, Ek);
}

/**
 * @brief Tính toán sai số đồng hồ của vệ tinh bao gồm cả hiệu ứng tương đối tính
 * @param Ek Dị thường lệch tâm (đã tính được từ vòng lặp Kepler)
 */
double calculate_sat_clock_bias(const Ephemeris* eph, double t, double Ek) {
    double tk = t - eph->t_oc;
    if (tk > 302400) tk -= 604800;
    else if (tk < -302400) tk += 604800;

    // Phương trình đa thức đồng hồ từ Subframe 1
    double bias = eph->a_f0 + eph->a_f1 * tk + eph->a_f2 * tk * tk;

    // Hiệu ứng tương đối tính (Relativistic correction factor)
    // F = -2 * sqrt(MU) / c^2 = -4.442807633e-10
    double F = -4.442807633e-10;
    double rel_corr = F * eph->e * eph->sqrtA * sin(Ek);

    // Sai số đồng hồ tổng hợp (trừ đi T_GD trễ nhóm phần cứng)
    return bias + rel_corr - eph->T_GD;
}