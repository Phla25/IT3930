// Code Loop Discriminator - Ước lượng độ lệch thời gian (pha mã) giữa mã PRN đầu vào và bản sao mã cục bộ
// Kết quả được dùng để điều khiển Code NCO, giúp căn chỉnh lại tốc độ sinh mã PRN
#ifndef CODE_LOOP_DISCRIMINATOR_H
#define CODE_LOOP_DISCRIMINATOR_H

#include "IntegrateAndDump.h"

// Tính toán sai số pha mã (ví dụ dùng thuật toán Non-coherent Early minus Late Power)
float CodeLoopDiscriminator_Evaluate(const CorrelatorOutputs *outputs);

#endif