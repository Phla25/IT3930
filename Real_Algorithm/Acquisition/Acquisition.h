// Acquisition algorithms for GPS signal acquisition
#ifndef ACQUISITION_H
#define ACQUISITION_H
#define PI 3.14159265358979323846f
#define DOPPLER_MAX 10000.0f
#define DOPPLER_STEP 500.0f
// Cấu trúc số phức 
typedef struct {
    float real; // Phần thực
    float imag; // Phần ảo
} Complex;

// Cấu trúc kết quả tìm kiếm
typedef struct {
    int is_acquired; // 1 nếu tín hiệu được tìm thấy, 0 nếu không
    float best_doppler; // Doppler shift tốt nhất
    int best_code_phase_index; // Code phase tốt nhất
    float max_correlation; // Giá trị tương quan tối đa
} AcquisitionResult;

#endif // ACQUISITION_H