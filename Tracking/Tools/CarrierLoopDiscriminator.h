// Carrier Loop Discriminator - Đo sai số pha giữa sóng mang đầu vào và bản sao sóng mang cục bộ
#ifndef CARRIER_LOOP_DISCRIMINATOR_H
#define CARRIER_LOOP_DISCRIMINATOR_H

#include "IntegrateAndDump.h"

// Tính toán sai số pha sóng mang từ 2 nhánh cực đại I_P và Q_P (ví dụ dùng hàm atan2 hoặc I_P * Q_P)
float CarrierLoopDiscriminator_Evaluate(const CorrelatorOutputs *outputs);

#endif