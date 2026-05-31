import numpy as np
import matplotlib.pyplot as plt
import os

# =====================================================================
# CẤU HÌNH
# =====================================================================
N_FFT            = 131072
NUMBER_OF_SAMPLES = 38192   # fs=38.192MHz * 1ms
DOPPLER_STEP     = 500      # Hz
DOPPLER_MAX      = 10000    # Hz
DISPLAY_HALF     = 2000     # Hiển thị ±2000 samples quanh đỉnh trong 3D

# =====================================================================

def plot_acquisition_3d():
    folder = os.path.join("Real_Algorithm", "Result_Acquisition")
    if not os.path.isdir(folder):
        print(f"[LỖI] Không tìm thấy thư mục: {folder}")
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
        print(f"Tổng samples đọc được : {len(data)}")
        print(f"Doppler bins thực tế  : {num_doppler}  (expected {expected_bins})")
        print(f"Dư (không dùng)       : {leftover} samples")

        if num_doppler == 0:
            print("  ⚠ File rỗng hoặc kích thước sai, bỏ qua.")
            continue

        matrix = data.reshape((num_doppler, N_FFT))

        # --- Lấy phần tín hiệu thực (bỏ zero-padding) ---
        mag = np.abs(matrix[:, :NUMBER_OF_SAMPLES])   # shape: (num_doppler, NUMBER_OF_SAMPLES)

        # --- Tìm đỉnh toàn cục ---
        flat_idx         = np.argmax(mag)
        peak_dopp_bin    = flat_idx // NUMBER_OF_SAMPLES
        peak_phase       = flat_idx %  NUMBER_OF_SAMPLES
        peak_dopp_hz     = -DOPPLER_MAX + peak_dopp_bin * DOPPLER_STEP
        peak_val         = mag[peak_dopp_bin, peak_phase]

        print(f"Đỉnh cao nhất        : phase={peak_phase} samples | Doppler bin={peak_dopp_bin} ({peak_dopp_hz:+.0f} Hz)")
        print(f"Giá trị đỉnh         : {peak_val:.2f}")

        # Normalize
        mag_norm = mag / peak_val

        doppler_axis = np.linspace(-DOPPLER_MAX, DOPPLER_MAX, num_doppler)

        # =============================================================
        # HÌNH 1 — 3D surface chỉ vùng ±DISPLAY_HALF quanh đỉnh
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
        # HÌNH 2 — 2 slice 2D để chẩn đoán
        # =============================================================
        fig2, axes = plt.subplots(1, 2, figsize=(14, 4))
        fig2.suptitle(f"2D Slices — {fname}", fontsize=12)

        # ---- Slice 1: Dọc theo Code Phase tại Doppler bin tốt nhất ----
        ax1 = axes[0]
        slice_code = mag_norm[peak_dopp_bin, :]        # toàn bộ code axis
        ax1.plot(np.arange(NUMBER_OF_SAMPLES), slice_code, lw=0.8, color='steelblue')
        ax1.axvline(x=peak_phase, color='red',    lw=1.5, ls='--', label=f'Peak = {peak_phase}')
        # Vẽ thêm điểm mirror để kiểm tra bug
        mirror_phase = (NUMBER_OF_SAMPLES - peak_phase) % NUMBER_OF_SAMPLES
        ax1.axvline(x=mirror_phase, color='orange', lw=1,   ls=':', label=f'Mirror = {mirror_phase}')
        ax1.set_title(f"Code Phase slice @ Doppler = {peak_dopp_hz:+.0f} Hz")
        ax1.set_xlabel("Code Phase (samples)")
        ax1.set_ylabel("Normalized Correlation")
        ax1.legend(fontsize=9)
        ax1.set_xlim(0, NUMBER_OF_SAMPLES)

        # ---- Slice 2: Dọc theo Doppler tại code phase tốt nhất ----
        ax2 = axes[1]
        slice_dopp = mag_norm[:, peak_phase]
        ax2.plot(doppler_axis, slice_dopp, lw=1.2, color='darkorange', marker='o', ms=3)
        ax2.axvline(x=peak_dopp_hz, color='red', lw=1.5, ls='--', label=f'{peak_dopp_hz:+.0f} Hz')
        ax2.set_title(f"Doppler slice @ Code Phase = {peak_phase} samples")
        ax2.set_xlabel("Doppler (Hz)")
        ax2.set_ylabel("Normalized Correlation")
        ax2.legend(fontsize=9)

        plt.tight_layout()
        plt.show()

        # =============================================================
        # CHẨN ĐOÁN: Có phải bug mirror không?
        # =============================================================
        val_at_mirror = mag_norm[peak_dopp_bin, mirror_phase]
        ratio_mirror  = val_at_mirror / 1.0   # so với đỉnh (đã normalize)
        print(f"\n[Chẩn đoán mirror]")
        print(f"  Giá trị tại mirror phase ({mirror_phase}): {val_at_mirror:.3f}")
        if val_at_mirror > 0.5:
            print("  ⚠ Có thể là BUG MIRROR — đỉnh thứ 2 tại vị trí đối xứng rất cao!")
            print("  → Kiểm tra lại phép tính actual_phase trong C.")
        else:
            print("  ✅ Không phải mirror — đỉnh thứ 2 là tín hiệu thật (multipath hoặc nhiễu).")

if __name__ == "__main__":
    plot_acquisition_3d()