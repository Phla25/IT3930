// Fine Frequency Search - implement trước khi đi vào tracking
#ifndef FINE_FREQUENCY_SEARCH_H
#define FINE_FREQUENCY_SEARCH_H
#include "Acquisition.h"
// Hàm thực hiện quá trình tìm kiếm tần số Doppler chính xác sau khi đã có kết quả Acquisition
// Cập nhật tham số để xử lý tín hiệu phức thay vì tín hiệu thực
AcquisitionResult performFineFrequencySearch(const Complex* signal_in, const float* local_prn, 
                                             int number_of_samples, float f_sampling, float if_f,
                                             AcquisitionResult coarse_result); 
#endif // FINE_FREQUENCY_SEARCH_H