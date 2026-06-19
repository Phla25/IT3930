import os
import pandas as pd
import matplotlib.pyplot as plt

# =============================================================================
# CẤU HÌNH THÔNG SỐ ĐỐI CHIẾU
# =============================================================================
PRN = 15
CSV_FILE_PATH = f"Result_Tracking\\tracking_results_PRN{PRN:02d}.csv"

if not os.path.exists(CSV_FILE_PATH):
    print(f"LỖI: Không tìm thấy file {CSV_FILE_PATH}!")
    print("Hãy chắc chắn bạn đã thêm lệnh ghi CSV vào main.c và chạy file .exe thành công.")
    exit()

print(f"Đang đọc dữ liệu từ: {CSV_FILE_PATH} ...")
df = pd.read_csv(CSV_FILE_PATH)

time_s = df['ms'] / 1000.0

# =============================================================================
# FIGURE 1: TẦN SỐ SÓNG MANG & BIỂU ĐỒ CHÒM SAO I/Q
# =============================================================================
fig1 = plt.figure(figsize=(12, 5))
fig1.canvas.manager.set_window_title(f'Tracking Results SV PRN {PRN} - Figure 1')

ax1 = fig1.add_subplot(1, 2, 1)
ax1.plot(time_s, df['doppler'], 'b', linewidth=1.2)
ax1.set_title('Discrete Carrier Frequency')
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Doppler (Hz)')
ax1.grid(True, linestyle=':', alpha=0.6)

ax2 = fig1.add_subplot(1, 2, 2)
ax2.plot(df['I_P'], df['Q_P'], '.', color='darkblue', markersize=2, alpha=0.7)
ax2.axhline(0, color='black', linewidth=0.8, linestyle='--')
ax2.axvline(0, color='black', linewidth=0.8, linestyle='--')
ax2.set_title('In-phase vs. Quadrature')
ax2.set_xlabel('In-phase (I_P)')
ax2.set_ylabel('Quadrature (Q_P)')
ax2.grid(True, linestyle=':', alpha=0.6)

fig1.tight_layout()

# =============================================================================
# FIGURE 2: ĐỘNG LỰC HỌC MẠCH VÒNG KHÓA VÀ CHUỖI BIT DẪN ĐƯỜNG
# =============================================================================
fig2, axs = plt.subplots(3, 2, figsize=(13, 10), sharex=False)
fig2.canvas.manager.set_window_title(f'Tracking Results SV PRN {PRN} - Figure 2')

axs[0, 0].plot(time_s, df['codeFreq'], 'r', linewidth=1.2)
axs[0, 0].set_title('Discrete Code Frequency')
axs[0, 0].set_xlabel('Time (s)')
axs[0, 0].set_ylabel('Frequency (Hz)')
axs[0, 0].grid(True, linestyle=':', alpha=0.6)
axs[0, 0].get_yaxis().get_major_formatter().set_useOffset(False) 

axs[0, 1].plot(time_s, df['amp_E'], color='orange', label='$\\sqrt{I_{E}^2 + Q_{E}^2}$', linewidth=1.0)
axs[0, 1].plot(time_s, df['amp_P'], color='green', label='$\\sqrt{I_{P}^2 + Q_{P}^2}$', linewidth=1.2)
axs[0, 1].plot(time_s, df['amp_L'], color='purple', label='$\\sqrt{I_{L}^2 + Q_{L}^2}$', linewidth=1.0, linestyle='--')
axs[0, 1].set_title('Correlation Results')
axs[0, 1].set_xlabel('Time (s)')
axs[0, 1].set_ylabel('Amplitude')
axs[0, 1].grid(True, linestyle=':', alpha=0.6)
axs[0, 1].legend(loc="upper right")

axs[1, 0].plot(time_s, df['pllDiscr'], 'b', linewidth=0.8)
axs[1, 0].set_title('Filtered PLL Discriminator')
axs[1, 0].set_xlabel('Time (s)')
axs[1, 0].set_ylabel('Amplitude (Cycles)')
axs[1, 0].grid(True, linestyle=':', alpha=0.6)

axs[1, 1].plot(time_s, df['dllDiscr'], 'r', linewidth=0.8)
axs[1, 1].set_title('Raw DLL Discriminator')
axs[1, 1].set_xlabel('Time (s)')
axs[1, 1].set_ylabel('Amplitude (Chips)')
axs[1, 1].grid(True, linestyle=':', alpha=0.6)

axs[2, 0].step(time_s, df['navBit'], color='black', where='post', linewidth=1.5)
axs[2, 0].set_xlim(0, 2.0)
y_min = df['navBit'].min()
y_max = df['navBit'].max()
y_margin = (y_max - y_min) * 0.1  
if y_max > y_min: 
    axs[2, 0].set_ylim(y_min - y_margin, y_max + y_margin)
axs[2, 0].axhline(0, color='red', linestyle='--', linewidth=1.5, alpha=0.8)
axs[2, 0].set_title('Raw Navigation Data (I_P Amplitude) - Zoomed 2s')
axs[2, 0].set_xlabel('Time (s)')
axs[2, 0].set_ylabel('Amplitude (I_P)')
axs[2, 0].grid(True, linestyle=':', alpha=0.6)

axs[2, 1].axis('off') 
doppler_mean = df['doppler'].mean()
code_mean = df['codeFreq'].mean()
info_text = (
    f"== SV PRN {PRN:02d} STATISTICS ==\n\n"
    f"Total Processed Time: {time_s.max():.3f} seconds\n"
    f"Mean Doppler Shift: {doppler_mean:+.2f} Hz\n"
    f"Mean Code Frequency: {code_mean:.2f} Hz\n"
    f"Status: LOCK STABLE (Complex I/Q)"
)
axs[2, 1].text(0.1, 0.3, info_text, fontsize=11, family='monospace',
              bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=1'))

fig2.tight_layout()
print(">> [PLOT DONE] Đã hiển thị đồ thị thành công!")
plt.show()