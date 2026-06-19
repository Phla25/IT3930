// Preamble.h - Định nghĩa mã Preamble và trạng thái máy (State Machine) cho việc đồng bộ khung dữ liệu GPS
#ifndef PREAMBLE_H
#define PREAMBLE_H
#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>

// Định nghĩa mã Preamble và dạng đảo pha của nó
#define PREAMBLE_NORMAL 0x8B // 10001011 trong hệ nhị phân
#define PREAMBLE_INVERT 0x74 // 01110100 (đảo ngược của 0x8B do lệch pha Costas)
// Trạng thái máy (State Machine)
typedef enum {
    SEARCHING_PREAMBLE,
    VERIFYING_PARITY,
    DECODING_SUBFRAME
} FrameSyncState;
void process_new_nav_bit(int new_bit);
int check_parity(const uint8_t* word_bits);
#endif // PREAMBLE_H