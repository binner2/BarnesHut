# Gerçek Verilerle Yapılacaklar

Durum (5 Kasım 2026): Bütün ölçümler sentetik sahnede alındı. Bu belge, gerçek
derinlik kamerası verisine geçiş için yapılacak işleri, sırasıyla ve kontrol
listesi hâlinde verir. Her madde için "çıktı" satırı, işin bittiğini kanıtlayan
dosyayı tanımlar.

## 0. Ön koşullar (1 gün)

- [ ] `libpng-dev`, `cmake`, `g++ >= 11`; derleme: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j && ctest --test-dir build`
- [ ] Python ortamı: `numpy pillow h5py open3d` (Open3D taban çizgileri için; sürümü not et)
- [ ] Zaman ölçümü için sabit donanım: CPU modeli, çekirdek sayısı, `OMP_NUM_THREADS`, governor `performance`; her CSV'ye makine kimliği eklenir
- [ ] Disk: TUM dizileri 0,5–2 GB/dizi, NYU v2 2,8 GB, ScanNet örnek taraması 1–3 GB

## 1. Veri kümeleri

### 1.1 TUM RGB-D (eksen A ve D; akış/çevrimiçi senaryo)

Kaynak: Sturm vd. 2012, doi:10.1109/IROS.2012.6385773. Biçim sayfası
(https://cvg.cit.tum.de/data/datasets/rgbd-dataset/file_formats, erişim 4 Kas 2026):
derinlik 640×480 16-bit PNG, **ölçek 5000** (ham/5000 = metre), 0 = eksik;
önerilen iç parametreler (ROS varsayılanı) fx = fy = 525,0, cx = 319,5, cy = 239,5.

İndirme kalıbı (indirme sayfası, erişim 5 Kas 2026):
`https://cvg.cit.tum.de/rgbd/dataset/freiburg{N}/rgbd_dataset_freiburg{N}_{dizi}.tgz`

| Dizi | Neden |
|------|-------|
| `freiburg1_xyz` | Sayfanın ilk deneyler için önerdiği dizi; küçük hareket, masa sahnesi |
| `freiburg1_desk`, `freiburg1_desk2` | Nesne yoğun masa; bölütleme (eksen C) |
| `freiburg2_xyz` | Daha geniş ofis, farklı kamera |
| `freiburg3_long_office_household` | Uzun dizi; akış/kare süresi dağılımı (eksen D) |
| `freiburg3_structure_notexture_far` | Düzlemsel, dokusuz; WLOP ve normal kestirimi (eksen B) |

Adımlar:
- [ ] Dizileri indir, `data/tum/<dizi>/depth/*.png` yapısını koru
- [ ] Her kare için: `pointcloud_tool depth2cloud depth/<ts>.png --fx 525 --fy 525 --cx 319.5 --cy 239.5 --scale 5000 -o clouds/<ts>.ply` (`tools/tum_depth_to_ply.sh`)
- [ ] fr1 ve fr2 için kalibre iç parametreler sayfada verilmiştir; "ROS varsayılanı" ile kalibre değerler arasındaki farkın KDE/WLOP sonuçlarına etkisini bir kez ölç ve rapora yaz
- [ ] Çıktı: `experiments/results/<tarih>/tum_<dizi>_frames.csv` (kare, nokta sayısı, inşa süresi, her algoritma süresi)

### 1.2 NYU Depth v2 (eksen C; etiketli bölütleme)

Kaynak: Silberman vd. 2012, doi:10.1007/978-3-642-33715-4_54. Etiketli küme:
`http://horatio.cs.nyu.edu/mit/silberman/nyu_depth_v2/nyu_depth_v2_labeled.mat`
(~2,8 GB; sayfa erişimi 5 Kas 2026). Sayfaya göre `depths` değişkeni H×W×N ve
**metre** birimindedir; dosya MATLAB v7.3 (HDF5) olduğundan `h5py` ile okunur.

- [ ] `tools/nyu_mat_to_depth_png.py` ile `depths` → 16-bit PNG (ölçek 1000, mm) ve `labels` → 16-bit PNG etiket haritası; iç parametreler NYU araç kutusundaki `camera_params.m` değerlerinden alınır (**bu değerler bu oturumda doğrulanmadı; araç kutusundan okunup buraya yazılacak**)
- [ ] Etiket haritası ile nokta bulutu etiketleri eşle (piksel → nokta indeksi; `depth_to_point_cloud` `pixel_of_point` çıktısını verir)
- [ ] Çıktı: `data/nyu/{depth,label}/NNNN.png`, `experiments/results/<tarih>/nyu_segmentation.csv` (kare, yöntem, mIoU, düzlem sayısı, süre)

### 1.3 ScanNet (eksen B; yeniden yapılandırılmış yüzey referansı)

Kaynak: Dai vd. 2017, doi:10.1109/CVPR.2017.261. Erişim, kullanım koşulları
formu ile GitHub `ScanNet/ScanNet` deposu üzerinden talep edilir; indirme
betiği ve biçimler orada belgelenmiştir (bu oturumda erişilmedi).

- [ ] 5–10 tarama için `.sens` derinlik karelerini ve `_vh_clean_2.ply` yüzeyini al
- [ ] Yüzey ağı, WLOP ve seyreltme çıktıları için Chamfer/Hausdorff referansı olur
- [ ] Çıktı: `experiments/results/<tarih>/scannet_resampling.csv`

### 1.4 Kendi kaydı (genellenebilirlik, canlı demo)

- [ ] RealSense D435/D455: `tools/realsense_capture.py` (pyrealsense2; **bu ortamda test edilmedi**) 16-bit PNG + iç parametre dosyası yazar; derinlik ölçeği cihazdan okunur (D4xx için tipik 0,001 m/birim ama cihazdan sorgulanmalı)
- [ ] Azure Kinect: `k4a` SDK ile aynı biçim; iç parametreler kalibrasyon blob'undan
- [ ] En az üç sahne: düzlemsel oda köşesi (WLOP, normal), nesneli masa (bölütleme), uzak koridor (gürültü modeli \(\sigma_z = c z^2\) kalibrasyonu)
- [ ] Gürültü kalibrasyonu: sabit kameradan 100 kare düz duvar; piksel başına derinlik standart sapması ↔ uzaklık eğrisine karesel uydurma; `c` değeri `synthetic` seçeneği `--noise` ile aynı birimde raporlanır

## 2. Taban çizgileri (karşılaştırma yöntemleri)

| Eksen | Taban | Nasıl |
|------|-------|-------|
| A | Yarıçap kesmeli kaba kuvvet | `pointcloud_tool density --theta 0 --cutoff 4` (mevcut) |
| A | Open3D `compute_nearest_neighbor_distance`, KD-tree yarıçap KDE | Python betiği; aynı h |
| A | IFGT / FIGTree | C++ kütüphanesi (Morariu vd.); kurulum ve sürüm kaydedilecek |
| B | Voxel centroid, FPS, Poisson-disk | Open3D `voxel_down_sample`, `farthest_point_down_sample`, `sample_points_poisson_disk` (ağ gerekir) |
| B | CGAL `wlop_simplify_and_regularize_point_set` | Aynı h, μ, yineleme |
| C | PCL RANSAC düzlem bölütleme, bölge büyütme | `pcl_segmentation`; Open3D `segment_plane` ardışık |
| D | Bu deponun CPU sürümü | `--threads 1` ve `--threads N` |

## 3. Ölçüm protokolü

- Her ayar için 5 tekrar; medyan ve 5–95 persentil raporlanır (sentetik benchmark 2 tekrarın en iyisini kullanıyordu; gerçek veride medyana geç)
- Süre, octree inşası dahil ve hariç ayrı ayrı
- Doğruluk: KDE için kesin toplama göre göreli RMS ve en büyük hata; mean shift için etiket uyumu ve mIoU; WLOP için NN CoV, Chamfer, normal açısı hatası
- Her CSV satırı: makine kimliği, derleyici, commit SHA, veri kümesi, kare, parametreler
- Şekiller `experiments/plot_benchmarks.py` ile CSV'den üretilir; elle çizim yok

## 4. Gerçek veriye geçmeden önce koddaki bilinen eksikler

- [ ] Değişken bant genişliği \(h(z)\): derinlikle büyüyen h; `KernelDensityOptions` ve WLOP'a `bandwidth_per_metre` ekle
- [ ] Mean shift için konum–normal 6-D çekirdek (eksen C); düğümlerde normal ortalaması ve dağılımı
- [ ] WLOP çekim teriminde kesme yarıçapı h yerine 1,5h denenmeli; yinelemeler arası örnek ağacının yerinde güncellenmesi
- [ ] Mean shift yineleme üst sınırı yerine adım normuna göre erken durdurma (düzlemlerde yavaş yakınsama)
- [ ] PCD `binary_compressed` okuma (LZF) — PCL çıktılarının çoğu bu biçimdedir
- [ ] Organize (640×480 ızgaralı) bulutu koruyan mod; piksel komşuluğundan normal kestirimi tabanı
- [ ] `info` çıktısına yoğunluk–derinlik eğrisi (eq. \(\rho_A = f_x f_y \cos\alpha / z^2\)) doğrulaması

## 5. Teslim ölçütü

Bu plan, aşağıdaki üç dosya oluştuğunda tamamlanmış sayılır:
`experiments/results/<tarih>/tum_*.csv`, `nyu_segmentation.csv`,
`scannet_resampling.csv`; ve `docs/derinlik_kamerasi_barnes_hut_coworks.html`
raporunun "Sınırlılıklar" bölümünden "gerçek sensör verisi kullanılmamıştır"
cümlesi silindiğinde.
