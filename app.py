import streamlit as st
import geopandas as gpd
import rasterio
from rasterio.mask import mask
import numpy as np
from shapely.geometry import mapping
from rasterio.warp import reproject, Resampling
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import matplotlib.patheffects as patheffects
from pyproj import Transformer
import tempfile
import os
import plotly.express as px
import pandas as pd
import folium
from streamlit_folium import folium_static

# --- Konfigurasi Halaman Streamlit dan CSS Kustom ---
st.set_page_config(
    layout="wide",
    page_title="Aplikasi Kelas Kesesuaian Lahan Tanaman Manggis Di Kecamatan Jasinga",
    page_icon="🗺️",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .main-header {
        background: linear-gradient(90deg, #2E8B57, #228B22);
        padding: 2rem;
        border-radius: 10px;
        margin-bottom: 2rem;
        color: white;
        text-align: center;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
    }
    
    .main-header h1 {
        font-size: 2.5rem;
        font-weight: bold;
        margin-bottom: 0.5rem;
    }
    
    .main-header p {
        font-size: 1.2rem;
        opacity: 0.9;
    }
    
    .info-card {
        background: white;
        padding: 1.5rem;
        border-radius: 10px;
        border-left: 5px solid #2E8B57;
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
        margin: 1rem 0;
    }
    
    .legend-item {
        display: flex;
        align-items: center;
        margin: 0.5rem 0;
        padding: 0.5rem;
        background: #f8f9fa;
        border-radius: 5px;
    }
    
    .legend-color {
        width: 20px;
        height: 20px;
        border-radius: 3px;
        margin-right: 10px;
        border: 1px solid #ddd;
    }
    
    .parameter-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 1rem;
        border-radius: 8px;
        margin: 0.5rem 0;
    }
    
    .stButton > button {
        background: linear-gradient(90deg, #2E8B57, #228B22);
        color: white;
        border: none;
        border-radius: 25px;
        padding: 0.75rem 2rem;
        font-weight: bold;
        transition: all 0.3s ease;
    }
    
    .stButton > button:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 8px rgba(0, 0, 0, 0.2);
    }
    
    .sidebar .stSelectbox {
        margin-bottom: 1rem;
    }
    /* Style untuk pesan progres */
    .stProgress > div > div > div > div {
        background-color: #2E8B57;
    }
    .stAlert {
        border-radius: 8px;
    }
</style>
""", unsafe_allow_html=True)

# --- Header Aplikasi ---
st.markdown("""
<div class="main-header">
    <h1>🌱 Aplikasi Kelas Kesesuaian Lahan Tanaman Manggis Di Kecamatan Jasinga</h1>
    <p>Kecamatan Jasinga, Kabupaten Bogor, Jawa Barat</p>
</div>
""", unsafe_allow_html=True)

# --- Sidebar untuk Navigasi dan Kontrol ---
with st.sidebar:
    st.markdown("### 🛠️ Panel Kontrol")
    
    selected_layer = st.selectbox(
        "📊 Pilih Bagian Aplikasi:",
        ["Beranda", "Peta Kesesuaian Lahan", "Informasi Parameter", "Statistik"],
        index=0
    )
    
    st.markdown("---")
    
    st.markdown("### 📋 Parameter Analisis (Bobot)")
    parameters = {
        "🌡️ Suhu": "15%",
        "🌧️ Curah Hujan": "15%",
        "🧪 pH Tanah": "15%",
        "🏔️ Ketinggian": "15%",
        "📐 Kemiringan": "15%",
        "🌿 Tutupan Lahan": "15%",
        "🪨 Tekstur Tanah": "10%"
    }
    
    for param, weight in parameters.items():
        st.markdown(f"""
        <div class="parameter-card">
            <strong>{param}</strong><br>
            <small>Bobot: {weight}</small>
        </div>
        """, unsafe_allow_html=True)

# --- Fungsi Cache untuk Pemrosesan Data Spasial ---
@st.cache_data(show_spinner=False)
def generate_suitability_map():
    """
    Menghasilkan peta kesesuaian lahan dengan menggabungkan dan mengklasifikasikan
    berbagai layer raster berdasarkan bobot yang ditentukan.
    """
    shapefile_path = "Data/Shp_jasingan.shp"
    try:
        if not os.path.exists(shapefile_path):
            st.error(f"❌ File shapefile tidak ditemukan: {shapefile_path}. Pastikan 'Shp_jasingan.shp' dan file pendukungnya ada di folder 'Data/'.")
            return None, None, None, None, None, None
        
        gdf = gpd.read_file(shapefile_path)
        gdf_jasinga = gdf[gdf['WADMKC'].str.lower() == 'jasinga']
        if gdf_jasinga.empty:
            st.error("⚠️ Wilayah 'Jasinga' tidak ditemukan dalam shapefile 'Shp_jasingan.shp'. Pastikan nama kolom 'WADMKC' dan nilai 'jasinga' sudah benar.")
            return None, None, None, None, None, None
    except Exception as e:
        st.error(f"❌ Error saat memuat shapefile: {e}. Pastikan file shapefile lengkap dan tidak rusak.")
        return None, None, None, None, None, None

    raster_files = {
        "suhu": "Data/Skor_Suhu_Jasinga.tif",
        "curah_hujan": "Data/Skor Curah Hujan.tif",
        "ph_tanah": "Data/Skor_PH_soil_0_5cm_kec_jasinga.tif",
        "tekstur_tanah": "Data/Skor_Tekstur_Tanah_Jasinga.tif",
        "ketinggian": "Data/Skor Ketinggian.tif",
        "kemiringan": "Data/Skor_Kemiringan_Jasingan.tif",
        "tutupan_lahan": "Data/Skor_Tutupan_Lahan_Jasinga.tif"
    }

    weights = {
        "suhu": 0.15,
        "curah_hujan": 0.15,
        "ph_tanah": 0.15,
        "tekstur_tanah": 0.10,
        "ketinggian": 0.15,
        "kemiringan": 0.15,
        "tutupan_lahan": 0.15
    }

    for name, path in raster_files.items():
        if not os.path.exists(path):
            st.error(f"❌ File raster tidak ditemukan: {path}. Pastikan semua file .tif ada di folder 'Data/'.")
            return None, None, None, None, None, None

    ref_path = list(raster_files.values())[0]
    try:
        with rasterio.open(ref_path) as ref:
            ref_crs = ref.crs
            ref_transform = ref.transform
            ref_shape = (ref.height, ref.width)
            ref_meta = ref.meta.copy()
    except Exception as e:
        st.error(f"❌ Error saat membuka raster referensi ({ref_path}): {e}. Pastikan file ini tidak rusak.")
        return None, None, None, None, None, None

    skor_total = np.zeros((ref_shape[0], ref_shape[1]), dtype='float32')
    out_meta = ref_meta.copy()
    out_meta.update({"dtype": "float32", "count": 1, "nodata": np.nan})

    progress_bar = st.progress(0)
    status_text = st.empty()
    
    total_rasters = len(raster_files)
    processed_rasters = 0

    for name, path in raster_files.items():
        status_text.text(f"Memproses layer: {name.replace('_', ' ').title()}...")
        weight = weights[name]
        try:
            with rasterio.open(path) as src:
                gdf_masked = gdf_jasinga.to_crs(src.crs) if gdf_jasinga.crs != src.crs else gdf_jasinga
                shapes = [mapping(geom) for geom in gdf_masked.geometry]

                try:
                    masked, _ = mask(src, shapes, crop=True, filled=True)
                except ValueError:
                    st.warning(f"⚠️ Tidak dapat memproses layer {name}. Mungkin tidak ada tumpang tindih antara shapefile dan raster ini.")
                    continue

                masked = masked.astype('float32')
                nodata_value = src.nodata if src.nodata is not None else np.nan
                masked[masked == nodata_value] = np.nan

                resampled = np.empty((ref_shape[0], ref_shape[1]), dtype='float32')
                reproject(
                    source=masked,
                    destination=resampled,
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=ref_transform,
                    dst_crs=ref_crs,
                    resampling=Resampling.bilinear,
                    src_nodata=nodata_value,
                    dst_nodata=np.nan
                )
                
                skor_total = np.nansum([skor_total, resampled * weight], axis=0)

            processed_rasters += 1
            progress_bar.progress(processed_rasters / total_rasters)

        except Exception as e:
            st.error(f"❌ Error saat memproses layer '{name}': {e}")
            return None, None, None, None, None, None

    status_text.text("✅ Proses perhitungan selesai!")
    progress_bar.empty()

    if gdf_jasinga.crs != ref_crs:
        gdf_jasinga = gdf_jasinga.to_crs(ref_crs)
    shapes = [mapping(geom) for geom in gdf_jasinga.geometry]

    with tempfile.TemporaryDirectory() as tmpdir:
        temp_raster_path = os.path.join(tmpdir, "temp_skor.tif")
        with rasterio.open(temp_raster_path, 'w', **out_meta) as dst:
            dst.write(skor_total, 1)

        with rasterio.open(temp_raster_path) as src_final_mask:
            final_masked_array, final_masked_transform = mask(src_final_mask, shapes, crop=True, nodata=np.nan)
            
    skor_total_masked = final_masked_array[0]

    klasifikasi = np.full(skor_total_masked.shape, np.nan)
    klasifikasi[(skor_total_masked >= 2.9)] = 4
    klasifikasi[(skor_total_masked >= 2.5) & (skor_total_masked < 2.9)] = 3
    klasifikasi[(skor_total_masked >= 1.8) & (skor_total_masked < 2.5)] = 2
    klasifikasi[(skor_total_masked > 0) & (skor_total_masked < 1.8)] = 1
    klasifikasi[np.isnan(skor_total_masked)] = np.nan

    return klasifikasi, gdf_jasinga, ref_transform, ref_crs, skor_total_masked.shape, final_masked_transform

# --- Fungsi untuk Membuat Visualisasi Peta Interaktif (Folium) ---
def create_suitability_map_folium(klasifikasi_data, gdf_jasinga, raster_transform):
    colors_map = {
        1: '#FF0000',  # Red (N - Tidak Sesuai)
        2: '#FFFF00',  # Yellow (S3 - Sesuai Marginal)
        3: '#90EE90',  # Lightgreen (S2 - Cukup Sesuai)
        4: '#008000'   # Green (S1 - Sangat Sesuai)
    }

    rgb_data = np.zeros((*klasifikasi_data.shape, 4), dtype=np.uint8)
    for val, hex_color in colors_map.items():
        r, g, b = colors.hex2color(hex_color)
        mask = (klasifikasi_data == val)
        rgb_data[mask, 0] = int(r * 255)
        rgb_data[mask, 1] = int(g * 255)
        rgb_data[mask, 2] = int(b * 255)
        rgb_data[mask, 3] = 204
    rgb_data[np.isnan(klasifikasi_data), 3] = 0

    height, width = klasifikasi_data.shape
    xmin, ymax = raster_transform * (0, 0)
    xmax, ymin = raster_transform * (width, height)
    
    from_crs = gdf_jasinga.crs
    to_crs = "EPSG:4326"
    transformer = Transformer.from_crs(from_crs, to_crs, always_xy=True)
    
    lon1, lat1 = transformer.transform(xmin, ymin)
    lon2, lat2 = transformer.transform(xmax, ymax)
    bounds_latlon = [[min(lat1, lat2), min(lon1, lon2)], [max(lat1, lat2), max(lon1, lon2)]]

    center_lat = (bounds_latlon[0][0] + bounds_latlon[1][0]) / 2
    center_lon = (bounds_latlon[0][1] + bounds_latlon[1][1]) / 2

    m = folium.Map(location=[center_lat, center_lon], zoom_start=12, control_scale=True, tiles="OpenStreetMap")

    folium.TileLayer('OpenStreetMap', name='OpenStreetMap').add_to(m)
    folium.TileLayer('CartoDB positron', name='CartoDB Positron').add_to(m)
    folium.TileLayer('CartoDB dark_matter', name='CartoDB Dark Matter').add_to(m)
    folium.TileLayer(
        tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
        attr='Esri',
        name='Esri Satellite',
        overlay=False,
        control=True
    ).add_to(m)

    folium.raster_layers.ImageOverlay(
        image=rgb_data,
        bounds=bounds_latlon,
        opacity=0.8,
        name='Kesesuaian Lahan (Klasifikasi)'
    ).add_to(m)

    gdf_jasinga_wgs84 = gdf_jasinga.to_crs(to_crs)
    folium.GeoJson(
        gdf_jasinga_wgs84.__geo_interface__,
        name='Batas Administrasi Jasinga',
        style_function=lambda x: {
            'color': 'black',
            'weight': 2,
            'fillOpacity': 0
        }
    ).add_to(m)

    folium.LayerControl().add_to(m)

    legend_html = """
    <div style="
        position: fixed; 
        bottom: 50px; left: 50px; width: 180px; height: 180px; 
        border:2px solid grey; z-index:9999; font-size:14px;
        background-color:white;
        opacity:0.9;
        padding: 10px;
        border-radius: 8px;
        ">
        <h4 style="margin-top:0;">Kesesuaian Lahan</h4>
        <i style="background:#008000; width:18px; height:18px; float:left; margin-right:5px; border:1px solid #ccc;"></i> S1 (Sangat Sesuai)<br>
        <i style="background:#90EE90; width:18px; height:18px; float:left; margin-right:5px; border:1px solid #ccc;"></i> S2 (Cukup Sesuai)<br>
        <i style="background:#FFFF00; width:18px; height:18px; float:left; margin-right:5px; border:1px solid #ccc;"></i> S3 (Sesuai Marginal)<br>
        <i style="background:#FF0000; width:18px; height:18px; float:left; margin-right:5px; border:1px solid #ccc;"></i> N (Tidak Sesuai)<br>
    </div>
    """
    m.get_root().html.add_child(folium.Element(legend_html))

    return m

# --- Fungsi untuk Menghitung Statistik Kesesuaian Lahan ---
def calculate_statistics(klasifikasi):
    valid_data = klasifikasi[~np.isnan(klasifikasi)]
    if len(valid_data) == 0:
        return None
    
    stats = {}
    total_pixels = len(valid_data)
    
    class_definitions = {
        1: 'Tidak Sesuai (N)',
        2: 'Sesuai Marginal (S3)',
        3: 'Cukup Sesuai (S2)',
        4: 'Sangat Sesuai (S1)'
    }
    
    for class_val, class_name in class_definitions.items():
        count = np.sum(valid_data == class_val)
        percentage = (count / total_pixels) * 100 if total_pixels > 0 else 0
        stats[class_name] = {'count': count, 'percentage': percentage}
    
    return stats

# --- Area Konten Utama Aplikasi ---
if selected_layer == "Beranda":
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.markdown("""
        <div class="info-card">
            <h3>🎯 Tentang Aplikasi</h3>
            <p>Aplikasi ini dikembangkan untuk menyediakan analisis dan visualisasi kesesuaian lahan 
            di <b>Kecamatan Jasinga, Kabupaten Bogor</b>. Dengan memanfaatkan data geospasial dan 
            metode analisis multi-kriteria, aplikasi ini membantu mengidentifikasi area yang paling 
            sesuai untuk berbagai tujuan pemanfaatan lahan.</p>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown("""
        <div class="info-card">
            <h3>📊 Metodologi Analisis</h3>
            <p>Analisis kesesuaian lahan dilakukan menggunakan metode <b>Weighted Overlay (Pembobotan Bertingkat)</b>, 
            yaitu teknik analisis spasial yang menggabungkan beberapa lapisan data raster berdasarkan bobot kepentingan 
            masing-masing parameter. Proses ini melibatkan:</p>
            <ul>
                <li><b>Normalisasi Data:</b> Setiap parameter diubah menjadi skor standar (misal: 1-5) berdasarkan kesesuaiannya.</li>
                <li><b>Pembobotan:</b> Setiap parameter diberi bobot sesuai dengan tingkat pengaruhnya terhadap kesesuaian lahan.</li>
                <li><b>Overlay:</b> Skor parameter dikalikan dengan bobotnya dan dijumlahkan untuk mendapatkan skor total kesesuaian lahan.</li>
                <li><b>Klasifikasi:</b> Skor total diklasifikasikan ke dalam kelas-kelas kesesuaian (Sangat Sesuai, Cukup Sesuai, dst.).</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
        <div class="info-card">
            <h3>🏛️ Informasi Lokasi</h3>
            <ul>
                <li><strong>Kecamatan:</strong> Jasinga</li>
                <li><strong>Kabupaten:</strong> Bogor</li>
                <li><strong>Provinsi:</strong> Jawa Barat</li>
                <li><strong>Koordinat Umum:</strong> Sekitar 6°30' LS, 106°30' BT</li>
            </ul>
            <p>Kecamatan Jasinga merupakan wilayah yang memiliki karakteristik geografis beragam, 
            menjadikannya menarik untuk studi kesesuaian lahan.</p>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown("""
        <div class="info-card">
            <h3>🎨 Legenda Kelas Kesesuaian Lahan</h3>
        </div>
        """, unsafe_allow_html=True)
        
        legend_items = [
            ("#008000", "S1 - Sangat Sesuai"),
            ("#90EE90", "S2 - Cukup Sesuai"), 
            ("#FFFF00", "S3 - Sesuai Marginal"),
            ("#FF0000", "N - Tidak Sesuai")
        ]
        
        for color, label in legend_items:
            st.markdown(f"""
            <div class="legend-item">
                <div class="legend-color" style="background-color: {color};"></div>
                <span>{label}</span>
            </div>
            """, unsafe_allow_html=True)

elif selected_layer == "Peta Kesesuaian Lahan":
    st.markdown("### 🗺️ Peta Kelas Kesesuaian Lahan Kecamatan Jasinga (Interaktif)")
    st.info("💡 Gunakan fitur zoom (+/-), geser peta, dan opsi layer di kanan atas untuk interaksi.")
    
    placeholder_message = st.empty()
    placeholder_message.info("🔄 Sedang memproses data raster dan membuat peta. Ini mungkin memakan waktu beberapa menit. Mohon tunggu...")

    klasifikasi, gdf_jasinga, ref_transform, ref_crs, masked_shape, masked_transform = None, None, None, None, None, None
    try:
        klasifikasi, gdf_jasinga, ref_transform, ref_crs, masked_shape, masked_transform = generate_suitability_map()
        placeholder_message.empty()
    except Exception as e:
        placeholder_message.error(f"❌ Terjadi kesalahan saat memproses data peta: {e}. Silakan periksa file data Anda.")

    if klasifikasi is not None:
        try:
            m = create_suitability_map_folium(klasifikasi, gdf_jasinga, masked_transform)
            folium_static(m, width=1000, height=700)
            
            st.markdown("### 📈 Ringkasan Statistik Kesesuaian Lahan")
            stats = calculate_statistics(klasifikasi)
            
            if stats:
                cols = st.columns(4)
                class_order_for_metrics = ['Sangat Sesuai (S1)', 'Cukup Sesuai (S2)', 'Sesuai Marginal (S3)', 'Tidak Sesuai (N)']
                
                for i, class_name_full in enumerate(class_order_for_metrics):
                    if class_name_full in stats:
                        data = stats[class_name_full]
                        with cols[i]:
                            st.metric(
                                label=class_name_full.split(' (')[0],
                                value=f"{data['percentage']:.1f}%",
                                delta=f"Jumlah Piksel: {data['count']}"
                            )
            else:
                st.warning("Tidak ada data statistik yang tersedia untuk ditampilkan.")
            
            st.markdown("---")
            st.info("💡 Peta di atas bersifat interaktif. Anda bisa memperbesar, menggeser, dan memilih lapisan.")
                
        except Exception as e:
            st.error(f"❌ Gagal menampilkan peta interaktif: {e}. Coba periksa koneksi internet Anda atau restart aplikasi.")
            st.info("Pesan kesalahan ini bisa terjadi jika ada masalah dengan koordinat atau rendering Folium.")
    else:
        st.warning("Peta tidak dapat ditampilkan karena masalah pemrosesan data di atas.")

elif selected_layer == "Informasi Parameter":
    st.markdown("### 📋 Detail Parameter Analisis")
    st.markdown("<p>Bagian ini menjelaskan setiap parameter yang digunakan dalam analisis kesesuaian lahan, termasuk deskripsi, satuan, bobot yang diberikan, dan sumber datanya.</p>", unsafe_allow_html=True)
    
    param_details = {
        "🌡️ Suhu": {
            "deskripsi": "Suhu rata-rata tahunan yang mempengaruhi pertumbuhan dan metabolisme tanaman. Suhu optimal diperlukan untuk proses fotosintesis dan respirasi yang efisien.",
            "satuan": "°C",
            "bobot": "15%",
            "sumber": "Data klimatologi regional (misalnya dari citra satelit atau stasiun cuaca terdekat)."
        },
        "🌧️ Curah Hujan": {
            "deskripsi": "Jumlah curah hujan tahunan, penting untuk ketersediaan air bagi tanaman dan drainase. Curah hujan yang memadai dan terdistribusi merata sangat krusial.",
            "satuan": "mm/tahun", 
            "bobot": "15%",
            "sumber": "Data klimatologi regional (misalnya dari stasiun cuaca, TRMM, CHIRPS)."
        },
        "🧪 pH Tanah": {
            "deskripsi": "Tingkat keasaman atau kebasaan tanah yang sangat mempengaruhi ketersediaan nutrisi dan aktivitas mikroorganisme. Rentang pH yang sesuai penting untuk penyerapan unsur hara oleh tanaman.",
            "satuan": "pH unit",
            "bobot": "15%", 
            "sumber": "Analisis laboratorium tanah dari sampel lapangan, atau peta tanah nasional/regional."
        },
        "🏔️ Ketinggian": {
            "deskripsi": "Ketinggian dari permukaan laut, memengaruhi suhu, tekanan atmosfer, dan jenis vegetasi alami. Ketinggian tertentu mungkin lebih cocok untuk jenis tanaman spesifik.",
            "satuan": "meter dpl",
            "bobot": "15%",
            "sumber": "Digital Elevation Model (DEM) seperti SRTM, ASTER GDEM."
        },
        "📐 Kemiringan": {
            "deskripsi": "Derajat kemiringan permukaan lahan, berhubungan langsung dengan risiko erosi, drainase, dan kemudahan pengolahan lahan. Lahan yang terlalu miring cenderung rawan erosi.",
            "satuan": "derajat",
            "bobot": "15%",
            "sumber": "Analisis dari Digital Elevation Model (DEM)."
        },
        "🌿 Tutupan Lahan": {
            "deskripsi": "Jenis penutupan permukaan bumi (misalnya hutan, sawah, permukiman), mengindikasikan penggunaan lahan saat ini dan potensi pengembangannya serta implikasi ekologisnya.",
            "satuan": "kategori",
            "bobot": "15%",
            "sumber": "Citra satelit resolusi tinggi (misalnya Landsat, Sentinel) dan klasifikasi tutupan lahan."
        },
        "🪨 Tekstur Tanah": {
            "deskripsi": "Komposisi relatif partikel pasir, debu, dan liat dalam tanah, mempengaruhi kapasitas menahan air, aerasi, dan drainase. Tekstur tanah yang ideal mendukung pertumbuhan akar tanaman.",
            "satuan": "kategori (misal: Lempung, Pasir, Debu, Lempung Berpasir, dll.)",
            "bobot": "10%",
            "sumber": "Peta tanah detail, hasil survei dan analisis lapangan."
        }
    }
    
    for param, details in param_details.items():
        with st.expander(f"**{param}** - Bobot: {details['bobot']}"):
            col1, col2 = st.columns(2)
            with col1:
                st.write(f"**Deskripsi:** {details['deskripsi']}")
                st.write(f"**Satuan:** {details['satuan']}")
            with col2:
                st.write(f"**Bobot:** {details['bobot']}")
                st.write(f"**Sumber Data:** {details['sumber']}")

elif selected_layer == "Statistik":
    st.markdown("### 📊 Analisis Statistik Kesesuaian Lahan")
    st.markdown("<p>Visualisasi ini menampilkan distribusi spasial kelas kesesuaian lahan dalam bentuk grafik yang mudah dipahami.</p>", unsafe_allow_html=True)
    
    placeholder_message_stats = st.empty()
    placeholder_message_stats.info("🔄 Menghitung statistik kesesuaian lahan. Mohon tunggu...")

    klasifikasi = None
    result_stats = generate_suitability_map()
    if result_stats[0] is not None:
        klasifikasi = result_stats[0]
        placeholder_message_stats.empty()
    else:
        placeholder_message_stats.error("❌ Gagal mendapatkan data untuk statistik. Silakan cek bagian 'Peta Kesesuaian Lahan' untuk detail error.")

    if klasifikasi is not None:
        stats = calculate_statistics(klasifikasi)
        
        if stats:
            df_stats_list = []
            class_order_for_graphs = ['Tidak Sesuai (N)', 'Sesuai Marginal (S3)', 'Cukup Sesuai (S2)', 'Sangat Sesuai (S1)']
            display_names = {
                'Tidak Sesuai (N)': 'Tidak Sesuai',
                'Sesuai Marginal (S3)': 'Sesuai Marginal',
                'Cukup Sesuai (S2)': 'Cukup Sesuai',
                'Sangat Sesuai (S1)': 'Sangat Sesuai'
            }

            for class_full_name in class_order_for_graphs:
                if class_full_name in stats:
                    data = stats[class_full_name]
                    df_stats_list.append({
                        'Kelas': display_names[class_full_name], 
                        'Persentase': data['percentage']
                    })
            df_stats = pd.DataFrame(df_stats_list)
            
            df_stats['Kelas'] = pd.Categorical(df_stats['Kelas'], 
                                               categories=[display_names[k] for k in class_order_for_graphs], 
                                               ordered=True)
            df_stats = df_stats.sort_values('Kelas')

            plot_colors = ['#FF0000', '#FFFF00', '#90EE90', '#008000']
            color_map = {
                display_names['Tidak Sesuai (N)']: plot_colors[0],
                display_names['Sesuai Marginal (S3)']: plot_colors[1],
                display_names['Cukup Sesuai (S2)']: plot_colors[2],
                display_names['Sangat Sesuai (S1)']: plot_colors[3]
            }
            
            col1, col2 = st.columns(2)
            
            with col1:
                fig_pie = px.pie(df_stats, values='Persentase', names='Kelas', 
                                 title='Distribusi Area Berdasarkan Kelas Kesesuaian Lahan',
                                 color='Kelas',
                                 color_discrete_map=color_map,
                                 hole=0.3
                                 )
                fig_pie.update_traces(textinfo='percent+label', pull=[0.05 if df_stats['Kelas'].iloc[i] == display_names[max(stats, key=lambda k: stats[k]['percentage'])] else 0 for i in range(len(df_stats))])
                st.plotly_chart(fig_pie, use_container_width=True)
            
            with col2:
                fig_bar = px.bar(df_stats, x='Kelas', y='Persentase',
                                 title='Persentase Area untuk Setiap Kelas Kesesuaian',
                                 color='Kelas',
                                 color_discrete_map=color_map,
                                 labels={'Persentase': 'Persentase Area (%)', 'Kelas': 'Kelas Kesesuaian'}
                                 )
                fig_bar.update_yaxes(rangemode="tozero")
                fig_bar.update_traces(marker_line_width=1, marker_line_color='black')
                st.plotly_chart(fig_bar, use_container_width=True)
            
            st.markdown("### 📋 Ringkasan Data Statistik")
            st.dataframe(df_stats.style.format({'Persentase': '{:.2f}%'}), use_container_width=True)
            
            if not df_stats.empty:
                most_suitable_row = df_stats.loc[df_stats['Persentase'].idxmax()]
                st.success(f"💡 **Wawasan Utama:** Berdasarkan analisis, kelas **'{most_suitable_row['Kelas']}'** merupakan kategori kesesuaian lahan yang paling dominan di Kecamatan Jasinga, meliputi **{most_suitable_row['Persentase']:.1f}%** dari total area yang teranalisis. Hal ini menunjukkan potensi atau tantangan utama di wilayah tersebut.")
            else:
                st.info("Tidak ada data untuk analisis wawasan.")
            
        else:
            st.warning("Tidak ada statistik yang dapat dihitung dari data yang tersedia.")
    else:
        st.warning("Data klasifikasi tidak tersedia untuk menghitung statistik. Mohon pastikan pemrosesan peta berhasil.")

# --- Footer Aplikasi ---
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #666; padding: 1rem; font-size: 0.9em;">
    <p>🌱 <strong>Aplikasi Kelas Kesesuaian Lahan Tanaman Manggis Di Kecamatan Jasinga</strong> 🌱</p>
    <p><em>Dikembangkan dengan ❤️ menggunakan Streamlit, GeoPandas, Rasterio, Folium, dan Plotly.</em></p>
    <p><i>Data per [Juni 2025] | Untuk tujuan demonstrasi dan pendidikan.</i></p>
</div>
""", unsafe_allow_html=True)