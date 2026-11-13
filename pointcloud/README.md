# Derinlik Kamerası Nokta Bulutu Pipeline'ı (Barnes-Hut tabanlı)

Bu dizin, repodaki Barnes-Hut octree fikrini derinlik kamerasından (Kinect,
RealSense, Azure Kinect, ToF) alınan nokta bulutlarına uygular. Temel gözlem
şudur: N-body kuvvet toplamı ile kernel yoğunluk kestirimi, mean-shift ve
WLOP yeniden örnekleme aynı cebirsel biçime sahiptir:

    S(q) = sum_j w_j K(q, p_j)

Barnes-Hut yaklaşımı, `q`'dan yeterince uzak (`genişlik / uzaklık < theta`)
bir octree hücresinin tüm noktalarını hücrenin ağırlıklı merkeziyle temsil
eder; toplam O(N^2) yerine O(N log N) maliyete iner ve `theta` parametresi
doğruluk/hız dengesini kontrol eder. `theta = 0` aynı kod yolundan kesin
(brute-force) sonucu verir; bu, hata analizinin referansıdır.

## Bileşenler

| Dosya | İçerik |
|-------|--------|
| `point_cloud.h/.cpp` | `PointCloud` (SoA: nokta, normal, renk, skaler), pinhole `CameraIntrinsics`, derinlik görüntüsü → nokta bulutu, PLY/PCD/XYZ okuma-yazma, 16-bit PNG/PGM/raw derinlik okuma |
| `point_octree.h/.cpp` | İndeks tabanlı (lineer) octree: bitişik düğüm dizisi, her düğüm bir indeks aralığına sahip; `radius_search`, `knn_search`, `traverse_bh` |
| `algorithms.h/.cpp` | Voxel downsample, PCA normal kestirimi (kameraya yönlendirme), SOR, **BH-KDE**, **BH-mean-shift**, **BH-WLOP**, en yakın komşu düzgünlük metrikleri |
| `synthetic.h/.cpp` | Sentetik iç mekân sahnesi (zemin, iki duvar, küre, kutu) için ışın izlemeli derinlik görüntüsü; `sigma_z = c z^2` gürültü modeli |
| `pointcloud_tool.cpp` | Komut satırı aracı |
| `tests/test_pointcloud.cpp` | Kendi kendini doğrulayan testler (`ctest` ile) |

Simülasyon ağacı (`../tree.h`) ile farklar: akıllı işaretçili düğümler yerine
`std::vector<OctreeNode>` ve `int32` çocuk indeksleri; her düğüm permüte
edilmiş indeks dizisinde `[begin, end)` aralığına sahiptir; düğüm başına heap
tahsisi yoktur. Oktant kodlaması ve açılma ölçütü orijinal kodla aynıdır.

## Derleme

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build
```

CUDA/OpenGL görselleştirici artık varsayılan olarak kapalıdır
(`-DBH_ENABLE_VISUALIZATION=ON` ile açılır). libpng bulunursa 16-bit PNG
derinlik görüntüleri doğrudan okunur; bulunmazsa `tools/depth_png_to_pgm.py`
ile PGM veya PLY'ye çevirin. `make` ile de derlenebilir (`make test`).

## Kullanım

```bash
# 1. Derinlik görüntüsünü nokta bulutuna çevir (TUM RGB-D: scale 5000)
build/pointcloud/pointcloud_tool depth2cloud depth.png \
    --fx 525 --fy 525 --cx 319.5 --cy 239.5 --scale 5000 -o frame.ply

# 2. Ön işleme
build/pointcloud/pointcloud_tool sor frame.ply --k 20 --std 2.0 -o clean.ply
build/pointcloud/pointcloud_tool normals clean.ply --k 20 --viewpoint 0,0,0 -o normals.ply
build/pointcloud/pointcloud_tool downsample clean.ply --voxel 0.01 -o ds.ply

# 3. Barnes-Hut kernel toplamları
build/pointcloud/pointcloud_tool density ds.ply --h 0.03 --theta 0.5 --compare -o density.ply
build/pointcloud/pointcloud_tool meanshift ds.ply --h 0.08 --theta 0.5 --min-size 100 -o segments.ply
build/pointcloud/pointcloud_tool resample ds.ply --count 5000 --h 0.08 --theta 0.5 -o uniform.ply

# 4. Sentetik veri ve benchmark (gerçek kamera olmadan)
build/pointcloud/pointcloud_tool synthetic -o synth.ply --pgm synth.pgm --labels
build/pointcloud/pointcloud_tool bench --sizes 20000,50000,100000 --csv bench.csv

# 5. Gerçek kare üzerinde ölçüm: istatistik + derinlik–NN profili, kare benchmark'ı (5 tekrar medyanı)
build/pointcloud/pointcloud_tool info frame.ply --csv info.csv --tag 001 --depth-profile prof.csv --bin 0.25 --fx 525
build/pointcloud/pointcloud_tool benchframe frame.ply --csv frames.csv --tag 001 --repeats 5 \
    --hs 0.02,0.05 --thetas 0.5,0.7 --kappas 0,1 --ms-thetas 0,0.5,0.7,1.0 --wlop-thetas 0,0.5,0.7,1.0
```

`benchframe`, `bench` ile aynı CSV şemasını kullanır ve satırları dosyaya ekler; `experiments/run_tum.sh`
bunu bir TUM dizisinin seçilmiş kareleri için çalıştırır.

Çıktı PLY dosyalarında `scalar` alanı yoğunluk ya da küme etiketini taşır ve
renkler buna göre atanır; CloudCompare veya MeshLab ile doğrudan açılabilir.

İç parametre dosyası biçimi (`--intrinsics k.txt`):

```
# fx fy cx cy [width height [depth_scale]]
525.0 525.0 319.5 239.5 640 480 5000
```

## Algoritma notları

* **KDE**: `rho(q) = 1/(N (2 pi)^{3/2} h^3) sum_j exp(-|q-p_j|^2 / 2h^2)`.
  Uzak hücre `n` için katkı `W_n exp(-|q-c_n|^2 / 2h^2)` alınır. `--cutoff s`
  ile `s h`'den uzak hücreler tamamen atlanır (Gauss çekirdeği için ihmal
  edilebilir katkı).
* **Mean shift**: `y <- sum_j p_j K / sum_j K`; uzak hücreler `W_n c_n K(|y-c_n|)`
  ile temsil edilir. Yakınsayan modlar `h/2` yarıçapında birleştirilir.
* **WLOP** (Huang vd. 2009): çekim terimi giriş bulutu octree'si üzerinden,
  itme terimi örnek bulutu octree'si üzerinden; her ikisi de `traverse_bh`
  ile. `theta(r) = exp(-16 r^2/h^2)`, `eta(r) = -r`, `mu` itme ağırlığı.
* Barnes-Hut açılma ölçütü: `2 * half / |q - centroid| < theta`; `q`'yu
  içeren hücre asla yaklaşıklanmaz (öz terim her zaman kesindir).
