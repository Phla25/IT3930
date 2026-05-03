#include <stdio.h>
#include "Tools.h"
#include <stdlib.h>
#include <math.h>
// Hàm đảo bit (bit-reversal) cho FFT
void bitReversal(Complex* data, int n) {
    int j = 0;
    for (int i = 0; i < n; i++){
        if (i < j) {
            Complex temp = data[i];
            data[i] = data[j];
            data[j] = temp;
        }
        int k = n / 2;
        while (k > 0 && k <= j) {
            j -= k;
            k /= 2;
        }
        j += k;
    }
}

// Triển khai thuật toán Cooley-Tukey Radix-2 FFT
void customFFT(const Complex* input, Complex* output, int n) {

    for (int i = 0; i < n; i++) {
        output[i] = input[i];
    }
    // Sắp xếp lại vị trí các phần tử (Bit Reversal)
    bitReversal(output, n);

    // Tính toán các nút FFT
    for (int size = 2; size <= n; size *= 2){
        int halfSize = size / 2;
        float angle = -2.0f * PI / size; // Góc cho nút FFT

        for (int i = 0; i < n; i += size){ // Chia khối theo từng size
            for (int j = 0; j < halfSize; j++){     // Duyệt qua nửa đầu của khối
                int even_index = i + j;
                int odd_index = i + j + halfSize;   // Lý do vì sao cần xếp lại vị trí trước đó (bit reversal) 
                                                    // Đảm bảo rằng odd_index luôn nằm trong nửa sau của khối hiện tại

                // Tính Twiddle factor: e^(-j * 2 * PI *j / size)
                float phase = j * angle;
                Complex twiddle = {cosf(phase), sinf(phase)};   // Tạo ra số phức đại diện cho góc xoay
                Complex temp = {
                    output[odd_index].real * twiddle.real - output[odd_index].imag * twiddle.imag,
                    output[odd_index].real * twiddle.imag + output[odd_index].imag * twiddle.real
                };  // Nhân số phức odd_index với twiddle factor để xoay nó về đúng vị trí trong mặt phẳng phức

                // Cập nhật giá trị FFT
                // y_{k+n/2} = x_k - twiddle * x_{k+n/2}
                // nhánh dưới sẽ nhận giá trị đã được xoay (đã nhân với twiddle) và trừ đi giá trị ở nhánh trên cũ (even_index)
                output[odd_index].real = output[even_index].real - temp.real;
                output[odd_index].imag = output[even_index].imag - temp.imag;
                // y_k = x_k + twiddle * x_{k+n/2}
                // nhánh trên sẽ nhận giá trị đã được xoay (đã nhân với twiddle) và cộng thêm vào giá trị ở nhánh trên cũ (even_index)
                output[even_index].real += temp.real;
                output[even_index].imag += temp.imag;
            }
        }
    }
}

// Hàm IFFT (đảo ngược FFT)
void customIFFT(const Complex* input, Complex* output, int n) {
    for (int i = 0; i < n; i++){
        output[i] = input[i];
    }
    // Sắp xếp lại vị trí các phần tử (Bit Reversal)
    bitReversal(output, n);

    // Tính toán các nút IFFT
    for (int size = 2; size <= n; size *= 2){
        int halfSize = size / 2;
        float angle = 2.0f * PI / size; // Góc cho nút IFFT - góc ngược với FFT

        for (int i = 0; i < n; i += size){ // Chia khối theo từng size
            for (int j = 0; j < halfSize; j++){     // Duyệt qua nửa đầu của khối
                int even_index = i + j;
                int odd_index = i + j + halfSize;   // Lý do vì sao cần xếp lại vị trí trước đó (bit reversal) 
                                                    // Đảm bảo rằng odd_index luôn nằm trong nửa sau của khối hiện tại

                // Tính Twiddle factor: e^(-j * 2 * PI *j / size)
                float phase = j * angle;
                Complex twiddle = {cosf(phase), sinf(phase)};   // Tạo ra số phức đại diện cho góc xoay
                Complex temp = {
                    output[odd_index].real * twiddle.real - output[odd_index].imag * twiddle.imag,
                    output[odd_index].real * twiddle.imag + output[odd_index].imag * twiddle.real
                };  // Nhân số phức odd_index với twiddle factor để xoay nó về đúng vị trí trong mặt phẳng phức

                // Cập nhật giá trị IFFT
                // y_{k+n/2} = x_k - twiddle * x_{k+n/2}
                // nhánh dưới sẽ nhận giá trị đã được xoay (đã nhân với twiddle) và trừ đi giá trị ở nhánh trên cũ (even_index)
                output[odd_index].real = output[even_index].real - temp.real;
                output[odd_index].imag = output[even_index].imag - temp.imag;
                // y_k = x_k + twiddle * x_{k+n/2}
                // nhánh trên sẽ nhận giá trị đã được xoay (đã nhân với twiddle) và cộng thêm vào giá trị ở nhánh trên cũ (even_index)
                output[even_index].real += temp.real;
                output[even_index].imag += temp.imag;
            }
        }
    }
    // Chuẩn hóa kết quả IFFT bằng cách chia cho n
    for (int i = 0; i < n; i++){
        output[i].real /= n;
        output[i].imag /= n;
    }
}