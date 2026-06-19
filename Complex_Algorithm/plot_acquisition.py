import numpy as np
import matplotlib.pyplot as plt
import os

# =====================================================================
# CẤU HÌNH CHO TÍN HIỆU PHỨC (5 MHz, Zero-IF)
# =====================================================================
N_FFT             = 16384   # Kích thước FFT đã hạ xuống để tối ưu
NUMBER_OF_SAMPLES = 5000    # fs = 5.0 MHz * 1ms = 5000 mẫu
DOPPLER_STEP      = 500     # Hz
DOPPLER_MAX       = 10000   # Hz
DISPLAY_HALF      = 500     # Hiển thị ±500 mẫu quanh đỉnh để zoom cho rõ

# =====================================================================

def plot_acquisition_3d():
    folder = "Result_Acquisition"
    if not os.path.isdir(folder):
        print(f"[LỖI] Không tìm thấy thư mục: {folder}. Hãy chắc chắn bạn đã un-comment hàm save_acq_bin trong PCPSA.c")
        return

    files = [f for f in os.listdir(folder) if f.endswith(".bin")]
    if not files:
        print("[LỖI] Không tìm thấy file .bin nào!")
        return

    expected_bins = int(2 * DOPPLER_MAX / DOPPLER_STEP) + 1  # = 41

    for fname in files:
        filepath = os.path.join(folder, fname)
        data = np.fromfile(filepath, dtype=np.complex64)

        num_doppler = len(data) // N_FFT
        leftover    = len(data) %  N_FFT

        print(f"\n{'='*60}")
        print(f"File : {fname}")
        print(f"Tổng mẫu đọc được : {len(data)}")
        print(f"Doppler bins      : {num_doppler}  (expected {expected_bins})")

        if num_doppler == 0:
            print("  ⚠ File rỗng hoặc kích thước sai, bỏ qua.")
            continue

        matrix = data.reshape((num_doppler, N_FFT))

        # Lấy phần biên độ của số phức (Magnitude)
        mag = np.abs(matrix[:, :NUMBER_OF_SAMPLES])   

        # Tìm đỉnh toàn cục
        flat_idx         = np.argmax(mag)
        peak_dopp_bin    = flat_idx // NUMBER_OF_SAMPLES
        peak_phase       = flat_idx %  NUMBER_OF_SAMPLES
        peak_dopp_hz     = -DOPPLER_MAX + peak_dopp_bin * DOPPLER_STEP
        peak_val         = mag[peak_dopp_bin, peak_phase]

        print(f"Đỉnh cao nhất     : phase = {peak_phase} mẫu | Doppler = {peak_dopp_hz:+.0f} Hz")

        # Chuẩn hóa về 1.0
        mag_norm = mag / peak_val
        doppler_axis = np.linspace(-DOPPLER_MAX, DOPPLER_MAX, num_doppler)

        # =============================================================
        # HÌNH 1 — 3D surface
        # =============================================================
        start = max(0, peak_phase - DISPLAY_HALF)
        end   = min(NUMBER_OF_SAMPLES, peak_phase + DISPLAY_HALF)
        mag_crop = mag_norm[:, start:end]
        code_crop = np.arange(start, end)
        X, Y = np.meshgrid(code_crop, doppler_axis)

        fig = plt.figure(figsize=(12, 7))
        ax  = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X, Y, mag_crop, cmap='viridis', alpha=0.95, edgecolor='none')
        ax.set_title(f"3D Acquisition — {fname}\nPeak @ phase={peak_phase}, doppler={peak_dopp_hz:+.0f} Hz")
        ax.set_xlabel('Code Phase (samples)')
        ax.set_ylabel('Doppler (Hz)')
        ax.set_zlabel('Normalized Correlation')
        plt.tight_layout()
        plt.show()

        # =============================================================
        # HÌNH 2 — 2D Slices
        # =============================================================
        fig2, axes = plt.subplots(1, 2, figsize=(14, 4))
        fig2.suptitle(f"2D Slices — {fname}", fontsize=12)

        ax1 = axes[0]
        slice_code = mag_norm[peak_dopp_bin, :]        
        ax1.plot(np.arange(NUMBER_OF_SAMPLES), slice_code, lw=0.8, color='steelblue')
        ax1.axvline(x=peak_phase, color='red', lw=1.5, ls='--', label=f'Peak = {peak_phase}')
        ax1.set_title(f"Code Phase slice @ Doppler = {peak_dopp_hz:+.0f} Hz")
        ax1.set_xlabel("Code Phase (samples)")
        ax1.set_ylabel("Normalized Correlation")
        ax1.legend(fontsize=9)
        ax1.set_xlim(0, NUMBER_OF_SAMPLES)

        ax2 = axes[1]
        slice_dopp = mag_norm[:, peak_phase]
        ax2.plot(doppler_axis, slice_dopp, lw=1.2, color='darkorange', marker='o', ms=3)
        ax2.axvline(x=peak_dopp_hz, color='red', lw=1.5, ls='--', label=f'{peak_dopp_hz:+.0f} Hz')
        ax2.set_title(f"Doppler slice @ Code Phase = {peak_phase}")
        ax2.set_xlabel("Doppler (Hz)")
        ax2.set_ylabel("Normalized Correlation")
        ax2.legend(fontsize=9)

        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    plot_acquisition_3d()