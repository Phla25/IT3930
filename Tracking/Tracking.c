// Tracking - Quản lý toàn bộ quá trình bám bắt tín hiệu GPS từ đầu vào thô đến xuất ra dữ liệu định vị
#include "Tracking.h"
#include "./Tools/PRNCodeGenerator.h"
#include "./Tools/NCOCarrierGenerator.h"
#include "./Tools/IntegrateAndDump.h"
#include "./Tools/CodeLoopDiscriminator.h"
#include "./Tools/CarrierLoopDiscriminator.h"
#include "./Tools/CarrierLoopFilter.h"
#include "./Tools/CodeLoopFilter.h" 
#include <stdio.h>
#include <stdlib.h>
#include <math.h> 

#define PI 3.14159265358979323846f
#define CARRIER_LOOP_BANDWIDTH 25.0f // Băng thông vòng khóa sóng mang (Hz)
#define CARRIER_LOOP_ZETA 0.707f     // Hệ số cản cho bộ lọc vòng khóa sóng mang
#define CODE_LOOP_BANDWIDTH 2.0f     // Băng thông vòng mã DLL (hẹp để lọc nhiễu tốt hơn)
#define CODE_LOOP_ZETA 0.707f        // Hệ số cản cho bộ lọc vòng mã DLL
#define T_INTEGRATION 0.001f         // Thời gian tích phân (1ms)

void Tracking_Init(TrackingChannel *channel, int prn, float fs, float if_freq, float init_doppler) {
    channel->prn = prn;
    channel->sampling_freq = fs;
    channel->base_if_freq = if_freq;

    // Khởi tạo bộ sinh mã C/A cho vệ tinh 
    PRNCodeGenerator_Init(&channel->code_state, prn, fs);

    // Khởi tạo bộ dao động nội NCO sóng mang
    float initial_carrier_freq = if_freq + init_doppler;
    NCOCarrierGenerator_Init(&channel->carrier_state, initial_carrier_freq);

    // Khởi tạo bộ lọc vòng lặp sóng mang (PLL/Costas)
    CarrierLoopFilter_Init(&channel->carrier_filter, CARRIER_LOOP_BANDWIDTH, CARRIER_LOOP_ZETA, T_INTEGRATION);

    channel->carrier_filter.old_output = init_doppler;
    
    // Khởi tạo bộ lọc vòng lặp mã (DLL)
    CodeLoopFilter_Init(&channel->code_filter, CODE_LOOP_BANDWIDTH, CODE_LOOP_ZETA, T_INTEGRATION);
    
    // Khởi tạo các thông số giám sát
    channel->current_doppler = init_doppler;
    channel->current_code_error = 0.0f;
    channel->current_carrier_error = 0.0f;
    channel->rem_code_phase = 0.0f;
    channel->rem_carr_phase = 0.0f;
}
// Chỉ thực hiện 1 khối, sau đó sẽ cộng dồn kết quả 20 khối để quyết định 1 bit dữ liệu định vị
void Tracking_ProcessBlock(TrackingChannel *channel, const float *signal_block, int* sample_offset, 
                           CorrelatorOutputs *outputs, float *nav_data_bit)
{
    // Tính block size xử lý trong 1ms (biến đổi)
    float code_phase_step = channel->code_state.code_freq / channel->sampling_freq;
    int blksize = (int)ceilf((1023.0f - channel->rem_code_phase) / code_phase_step);

    float *early_code  = (float*)malloc(blksize * sizeof(float));
    float *prompt_code = (float*)malloc(blksize * sizeof(float));
    float *late_code   = (float*)malloc(blksize * sizeof(float));
    float *cos_carrier = (float*)malloc(blksize * sizeof(float));
    float *sin_carrier = (float*)malloc(blksize * sizeof(float));

    if (early_code == NULL || prompt_code == NULL || late_code == NULL || 
        cos_carrier == NULL || sin_carrier == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(early_code);
        free(prompt_code);
        free(late_code);
        free(cos_carrier);
        free(sin_carrier);
        return;
    }

    // 1. Sinh mã PRN và Sóng mang dựa trên trạng thái NCO của chu kỳ hiện tại
    PRNCodeGenerator_Generate(&channel->code_state, blksize, channel->sampling_freq, 0.5f, 
                               early_code, prompt_code, late_code);
                               
    NCOCarrierGenerator_Generate(&channel->carrier_state, blksize, channel->sampling_freq, 
                                  cos_carrier, sin_carrier);

    // 2. Tích phân và xả (Integrate and Dump) để lấy 6 giá trị I/Q 
    IntegrateAndDump_Process(signal_block, blksize, cos_carrier, sin_carrier, 
                              early_code, prompt_code, late_code, outputs); 
                    
    //ĐÓNG VÒNG LẶP SÓNG MANG (CARRIER TRACKING FEEDBACK)
    // Đo đạc sai số và đẩy qua bộ lọc
    channel->current_carrier_error = CarrierLoopDiscriminator_Evaluate(outputs);
    float doppler_correction = CarrierLoopFilter_Update(&channel->carrier_filter, channel->current_carrier_error);

    // Phản hồi: Cập nhật lại tần số Carrier NCO cho chu kỳ tiếp theo [1]
    // Vòng phản hồi âm 
    channel->current_doppler = doppler_correction;
    channel->carrier_state.carrier_freq = channel->base_if_freq - channel->current_doppler;

    // ĐÓNG VÒNG LẶP MÃ (CODE TRACKING FEEDBACK)
    // Đo đạc sai số pha mã và đẩy qua bộ lọc vòng mã (Code Loop Filter)
    channel->current_code_error = CodeLoopDiscriminator_Evaluate(outputs);
    float code_correction = CodeLoopFilter_Update(&channel->code_filter, channel->current_code_error);

    // Kỹ thuật Carrier Aiding: Sử dụng Doppler của sóng mang để hỗ trợ vòng mã
    // Tỷ lệ chuyển đổi = Tần số chip (1.023 MHz) / Tần số L1 (1575.42 MHz)
    float carrier_aiding = channel->current_doppler * (1023000.0f / 1575420000.0f);
    
    // Phản hồi: Cập nhật tốc độ sinh mã cho Code NCO chu kỳ tiếp theo [1]
    channel->code_state.code_freq = 1023000.0f - carrier_aiding + code_correction;

    // 3. Trích xuất dữ liệu định vị thô (Lưu ý: cần cộng dồn 20ms ở hàm bên ngoài để quyết định 1 bit)
    *nav_data_bit = outputs->I_P;

    // 4 Lưu phần dư pha mã và sóng mang chưa xử lý để cộng dồn vào chu kỳ tiếp theo
    channel->rem_code_phase = channel->rem_code_phase 
                            + blksize * code_phase_step - 1023.0f;
    // Cập nhật offset ra ngoài
    *sample_offset += blksize;

    free(early_code);
    free(prompt_code);
    free(late_code);
    free(cos_carrier);
    free(sin_carrier);
}