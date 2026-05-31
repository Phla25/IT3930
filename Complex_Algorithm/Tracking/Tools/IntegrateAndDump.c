// Integrate and Dump - Công đoạn nhân tín hiệu thô với sóng mang nội và mã PRN, sau đó tích phân (cộng dồn) trong 1ms để chuẩn bị cho các bộ Discriminator xử lý
#include "IntegrateAndDump.h"

void IntegrateAndDump_Process(const Complex *incoming_signal, int number_of_samples,
                              const float *cos_carrier, const float *sin_carrier,
                              const float *early_code, const float *prompt_code, const float *late_code,
                              CorrelatorOutputs *outputs)
{
    
    float accummulate_IE = 0.0f, accummulate_IP = 0.0f, accummulate_IL = 0.0f;
    float accummulate_QE = 0.0f, accummulate_QP = 0.0f, accummulate_QL = 0.0f;

    for (int i = 0; i < number_of_samples; i++) {
        // 0. Mẫu tín hiệu thô hiện tại
        float s_I = incoming_signal[i].real;
        float s_Q = incoming_signal[i].imag;

        // 1. Giai đoạn hạ tần xuống băng tần cơ sở (Carrier Wipe-Off)
        // Áp dụng phép nhân số phức tiêu chuẩn vì sin_carrier đã ngậm dấu âm từ NCO
        // Real = (s_I * cos) - (s_Q * sin)
        // Imag = (s_I * sin) + (s_Q * cos)
        float baseband_I = (s_I * cos_carrier[i]) - (s_Q * sin_carrier[i]); 
        float baseband_Q = (s_I * sin_carrier[i]) + (s_Q * cos_carrier[i]); 

        // 2. Giai đoạn tương quan mã (Code Correlation) & Tích phân (Integrate)
        // Nhân dòng I và Q cơ sở với 3 pha mã dịch lệch và cộng dồn liên tục
        // Nhánh Early (Sớm)
        accummulate_IE += baseband_I * early_code[i];
        accummulate_QE += baseband_Q * early_code[i];

        // Nhánh Prompt (Đúng pha)
        accummulate_IP += baseband_I * prompt_code[i];
        accummulate_QP += baseband_Q * prompt_code[i];

        // Nhánh Late (Trễ)
        accummulate_IL += baseband_I * late_code[i];
        accummulate_QL += baseband_Q * late_code[i];
    }

    // 3. Giai đoạn Xả (Dump)
    // Đẩy toàn bộ năng lượng đã tích lũy sau 1ms vào cấu trúc đầu ra để các bộ Discriminator xử lý
    outputs->I_E = accummulate_IE;
    outputs->I_P = accummulate_IP;
    outputs->I_L = accummulate_IL;
    
    outputs->Q_E = accummulate_QE;
    outputs->Q_P = accummulate_QP;
    outputs->Q_L = accummulate_QL;
}