// Preamble.c - Xử lý đồng bộ khung dữ liệu GPS bằng cách phát hiện Preamble và kiểm tra Parity
#include "Preamble.h"
FrameSyncState current_state = SEARCHING_PREAMBLE;

// Dùng mảng 62 bit làm cửa sổ trượt liên tục
uint8_t word_buffer[62] = {0}; 
int bits_collected = 0;

void process_new_nav_bit(int new_bit) {
    if (current_state == SEARCHING_PREAMBLE) {
        
        // 1. Trượt toàn bộ cửa sổ sang trái 1 bit (Giống băng chuyền)
        for (int i = 0; i < 61; i++) {
            word_buffer[i] = word_buffer[i+1];
        }
        word_buffer[61] = new_bit; // Nạp bit mới tinh vào cuối
        
        if (bits_collected < 62) {
            bits_collected++;
            return; // Đợi đủ 62 bit mới bắt đầu kiểm tra
        }

        // 2. Soi Preamble liên tục trên 8 bit (từ vị trí 2 đến 9)
        uint8_t current_8bit = 0;
        for (int i = 2; i <= 9; i++) {
            current_8bit = (current_8bit << 1) | word_buffer[i];
        }

        // 3. Nếu thấy bóng dáng Preamble, test Parity ngay lập tức!
        if (current_8bit == PREAMBLE_NORMAL || current_8bit == PREAMBLE_INVERT) {
            // Test liền lúc 2 Word. word_buffer[0..31] và word_buffer[30..61]
            if (check_parity(word_buffer) && check_parity(word_buffer + 30)) {
                printf("  -> [Frame Sync] PARITY MATCH 100%%! Chinh thuc khoa khung.\n");
                current_state = DECODING_SUBFRAME;
            }
        }
    }
}
/*
 * Hàm kiểm tra Parity cho 1 Word GPS (30 bit)
 * Đầu vào: Mảng word_bits chứa 32 bit (giá trị 0 hoặc 1)
 * - word_bits[0]      : D_29* (Bit Parity thứ 5 của Word TRƯỚC ĐÓ)
 * - word_bits[1]      : D_30* (Bit Parity thứ 6 của Word TRƯỚC ĐÓ)
 * - word_bits[2..25]  : D_1 đến D_24 (24 bit Data nhận được)
 * - word_bits[26..31] : D_25 đến D_30 (6 bit Parity nhận được)
 * Đầu ra:
 * 1 : Parity ĐÚNG, cực tính chuẩn (Không cần lật bit)
 * -1 : Parity ĐÚNG, nhưng cực tính bị NGƯỢC (Báo hệ thống phải lật toàn bộ bit)
 * 0 : Parity SAI (Lỗi dữ liệu hoặc bắt nhầm Preamble)
 */
int check_parity(const uint8_t *word_bits) {
    uint8_t D29_star = word_bits[0];
    uint8_t D30_star = word_bits[1];

    // 1. Tính d1 đến d24 (Data bits nguyên thủy)
    // Công thức trong ảnh: D_i = d_i XOR D_30* ==> d_i = D_i XOR D_30*
    // Mục đích: Khôi phục lại chuỗi bit gốc nếu bị lật pha
    uint8_t d[25]; // Mảng 1-indexed (d[1] đến d[24]) để code khớp y hệt ảnh
    for(int i = 1; i <= 24; i++) {
        d[i] = word_bits[1 + i] ^ D30_star; 
    }

    // 2. Tính toán 6 bit Parity theo bảng phương trình (Sử dụng toán tử ^)
    uint8_t D25 = D29_star ^ d[1] ^ d[2] ^ d[3] ^ d[5] ^ d[6] ^ d[10] ^ d[11] ^ d[12] ^ d[13] ^ d[14] ^ d[17] ^ d[18] ^ d[20] ^ d[23];
    uint8_t D26 = D30_star ^ d[2] ^ d[3] ^ d[4] ^ d[6] ^ d[7] ^ d[11] ^ d[12] ^ d[13] ^ d[14] ^ d[15] ^ d[18] ^ d[19] ^ d[21] ^ d[24];
    uint8_t D27 = D29_star ^ d[1] ^ d[3] ^ d[4] ^ d[5] ^ d[7] ^ d[8] ^ d[12] ^ d[13] ^ d[14] ^ d[15] ^ d[16] ^ d[19] ^ d[20] ^ d[22];
    uint8_t D28 = D30_star ^ d[2] ^ d[4] ^ d[5] ^ d[6] ^ d[8] ^ d[9] ^ d[13] ^ d[14] ^ d[15] ^ d[16] ^ d[17] ^ d[20] ^ d[21] ^ d[23];
    uint8_t D29 = D30_star ^ d[1] ^ d[3] ^ d[5] ^ d[6] ^ d[7] ^ d[9] ^ d[10] ^ d[14] ^ d[15] ^ d[16] ^ d[17] ^ d[18] ^ d[21] ^ d[22] ^ d[24];
    uint8_t D30 = D29_star ^ d[3] ^ d[5] ^ d[6] ^ d[8] ^ d[9] ^ d[10] ^ d[11] ^ d[13] ^ d[15] ^ d[19] ^ d[22] ^ d[23] ^ d[24];

    // 3. So sánh Parity tính được với Parity nhận được từ vệ tinh
    if (D25 == word_bits[26] && D26 == word_bits[27] && D27 == word_bits[28] &&
        D28 == word_bits[29] && D29 == word_bits[30] && D30 == word_bits[31]) {
        
        // Khớp hoàn hảo! Kiểm tra xem mạch vòng có đang bị kẹt ở pha 180 độ không.
        // Nếu D30_star == 1, nghĩa là toàn bộ bit đang bị lật ngược.
        if (D30_star == 1) {
            return -1; // Trả về -1 để báo "Hợp lệ, nhưng hãy lật bit đi"
        } else {
            return 1;  // Trả về 1 báo "Hợp lệ, giữ nguyên bit"
        }
    }
    // Không khớp Parity
    return 0; 
}