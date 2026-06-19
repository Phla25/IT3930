import os
import re
import folium
import matplotlib.pyplot as plt
import pandas as pd

# 🌟 BẢN VÁ ĐƯỜNG DẪN TUYỆT ĐỐI: Ép Python tìm file text đúng nơi file .py đang sống
script_dir = os.path.dirname(os.path.abspath(__file__)) # Lấy thư mục thật của file plot_map.py
log_filename = os.path.join(script_dir, "console_output.txt") # Ghép nối chính xác ra file text

# Định nghĩa thư mục đầu ra và tạo tự động bên cạnh file .py
output_dir = os.path.join(script_dir, "Result_Navigation")
os.makedirs(output_dir, exist_ok=True) 

# Kiểm tra lại một lần nữa
if not os.path.exists(log_filename):
    print(f"[LỖI] Đã quét tại đường dẫn thực tế:\n👉 '{log_filename}'\nNhưng vẫn không tìm thấy file! Bạn hãy kiểm tra lại xem tên file có bị gõ sai chính tả không nhé.")
    exit()

lats, lons = [], []

# REGEX bóc tách dựa trên cấu trúc [EPOCH...], chấp mọi lỗi font mã hóa chữ tiếng Việt trên Terminal Windows
pattern = r"\[EPOCH t0 \+ \d+ sec\]\s*->\s*[^:]+:\s*([0-9.-]+)[^|]*\|\s*[^:]+:\s*([0-9.-]+)"

with open(log_filename, "r", encoding="utf-8", errors="ignore") as f:
    for line in f:
        match = re.search(pattern, line)
        if match:
            lats.append(float(match.group(1)))
            lons.append(float(match.group(2)))

if not lats:
    print("[LỖI] Không thể bóc tách dữ liệu! Hãy kiểm tra lại nội dung file log.")
    exit()

# Đóng gói dữ liệu vào DataFrame
df = pd.DataFrame({"Latitude": lats, "Longitude": lons})
mean_lat = df["Latitude"].mean()
mean_lon = df["Longitude"].mean()

print(f"✅ Đã quét thành công Đám mây gồm {len(df)} Kỷ nguyên (Epoch).")
print(f"📍 Tâm đám mây hội tụ (Vị trí ước lượng): {mean_lat:.8f}, {mean_lon:.8f}")

# =============================================================================
# ĐỒ THỊ PHÂN TÁN HÌNH HỌC (MATPLOTLIB SCATTER CLOUD)
# =============================================================================
plt.figure(figsize=(7, 7))
plt.scatter(df["Longitude"], df["Latitude"], color="red", alpha=0.7, edgecolors="darkred", s=100, label="Epoch Solutions")
plt.scatter(mean_lon, mean_lat, color="blue", marker="X", s=250, label="Estimated Center (Mean)")

plt.title("GNSS Positioning Scatter Cloud (CEP Analysis)", fontsize=12, fontweight='bold')
plt.xlabel("Longitude (Degrees)", fontsize=10)
plt.ylabel("Latitude (Degrees)", fontsize=10)
plt.grid(True, linestyle="--", alpha=0.5)
plt.legend(loc="upper right")
plt.ticklabel_format(useOffset=False)

# Lưu đồ thị vào thư mục Result_Navigation
plot_png_path = os.path.join(output_dir, "gnss_scatter_cloud.png")
plt.savefig(plot_png_path, dpi=300, bbox_inches='tight')
print(f"🖼️  Đã lưu đồ thị phân tán: '{plot_png_path}'")

# =============================================================================
# ÉP LÊN BẢN ĐỒ VỆ TINH INTERACTIVE (FOLIUM WEB MAP)
# =============================================================================
mymap = folium.Map(location=[mean_lat, mean_lon], zoom_start=19, tiles="OpenStreetMap")

# Nhúng lớp Bản đồ nền Vệ tinh thực tế độ phân giải cao của Esri ArcGIS
folium.TileLayer(
    tiles="https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
    attr="EsriWorldImagery",
    name="Bản đồ Vệ tinh (Esri)",
    overlay=False,
    control=True
).add_to(mymap)

folium.LayerControl().add_to(mymap)

for idx, row in df.iterrows():
    folium.CircleMarker(
        location=[row["Latitude"], row["Longitude"]],
        radius=5,
        color="red",
        fill=True,
        fill_color="red",
        fill_opacity=0.5,
        popup=f"Epoch t0 + {idx:02d} sec"
    ).add_to(mymap)

folium.Marker(
    location=[mean_lat, mean_lon],
    popup=f"Tâm đám mây: {mean_lat:.6f}, {mean_lon:.6f}",
    icon=folium.Icon(color="blue", icon="info-sign")
).add_to(mymap)

# Lưu file HTML bản đồ vào thư mục Result_Navigation
output_html_path = os.path.join(output_dir, "gnss_scatter_cloud_map.html")
mymap.save(output_html_path)
print(f"🌎 Đã xuất bản đồ vệ tinh: '{output_html_path}'")

print("\n🎉 Toàn bộ kết quả định vị đã nằm gọn trong thư mục Result_Navigation!")
plt.show()