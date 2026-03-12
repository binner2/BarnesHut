# CUDA 13.2 Özelliklerinin Barnes-Hut N-Body Simülasyonuna Uygulanması

## Kapsamlı Analiz Raporu

**Proje:** Barnes-Hut N-Body Simülasyonu (C++20 / CUDA)
**Hedef CUDA Sürümü:** 13.2 (Mart 2026)
**Mevcut CUDA Kullanımı:** Yalnızca görselleştirme (rendering) — çekirdek algoritma CPU/OpenMP
**Tarih:** 2026-03-12

---

## İçindekiler

1. [Yönetici Özeti](#1-yönetici-özeti)
2. [Özellik 1: CUDA Tile (cuTile) Desteği](#2-özellik-1-cuda-tile-cutile-desteği)
3. [Özellik 2: CCCL 3.2 Modern C++ API'leri](#3-özellik-2-cccl-32-modern-c-apileri)
4. [Özellik 3: cub::DeviceTopK ve Segmented Reduction](#4-özellik-3-cubdevicetopk-ve-segmented-reduction)
5. [Özellik 4: Geliştirilmiş Bellek Yönetimi API'leri](#5-özellik-4-geliştirilmiş-bellek-yönetimi-apileri)
6. [Özellik 5: C++20 Uyumluluk İyileştirmeleri (nvcc)](#6-özellik-5-c20-uyumluluk-iyileştirmeleri-nvcc)
7. [Özellik 6: CUDA Graphs Geliştirmeleri](#7-özellik-6-cuda-graphs-geliştirmeleri)
8. [Özellik 7: Host Task Spin-Wait Modu](#8-özellik-7-host-task-spin-wait-modu)
9. [Özellik 8: Matematik Kütüphanesi Optimizasyonları](#9-özellik-8-matematik-kütüphanesi-optimizasyonları)
10. [Özellik 9: Gelişmiş Mimari Desteği (Blackwell)](#10-özellik-9-gelişmiş-mimari-desteği-blackwell)
11. [Özellik 10: Birleşik ARM Toolkit](#11-özellik-10-birleşik-arm-toolkit)
12. [Geçiş Stratejisi ve Öncelik Matrisi](#12-geçiş-stratejisi-ve-öncelik-matrisi)
13. [Karşılaştırmalı Performans Tahmin Tablosu](#13-karşılaştırmalı-performans-tahmin-tablosu)

---

## 1. Yönetici Özeti

Barnes-Hut N-Body simülasyonu şu anda **çekirdek hesaplama** (ağaç oluşturma, kuvvet hesabı, integrasyon) için yalnızca CPU ve OpenMP kullanmaktadır. CUDA yalnızca görselleştirme katmanında vertex hazırlama ve istatistik hesaplama amacıyla kullanılmaktadır.

CUDA 13.2, bu projeye aşağıdaki kritik geçiş fırsatlarını sunmaktadır:

| Öncelik | Özellik | Beklenen Etki | Zorluk |
|---------|---------|---------------|--------|
| **P0** | CCCL 3.2 Modern API | Temel GPU hızlandırma altyapısı | Orta |
| **P0** | cub Segmented Reduction | Kütle merkezi hesabında 10-66x hız | Düşük |
| **P1** | CUDA Graphs | Simülasyon döngüsü overhead'ini %40-60 azaltma | Orta |
| **P1** | Gelişmiş Bellek API'leri | Asenkron transfer optimizasyonu | Düşük |
| **P2** | CUDA Tile | Kuvvet hesabı kernel'inde veri yeniden kullanımı | Yüksek |
| **P2** | C++20 nvcc uyumu | Mevcut kodun doğrudan GPU derlemesi | Düşük |
| **P3** | Blackwell mimarisi | Yeni nesil donanım desteği | Düşük |
| **P3** | Matematik optimizasyonları | expm1f/erff hızlanması | Düşük |

---

## 2. Özellik 1: CUDA Tile (cuTile) Desteği

### Ne?
CUDA Tile, artık Compute Capability 8.x (Ampere, Ada) ve 10.x/12.x (Blackwell) mimarilerinde desteklenmektedir. Tile-tabanlı programlama, shared memory üzerinde yapılandırılmış veri erişim kalıpları sağlar.

### Neden Bu Projeye Uygulanmalı?

Barnes-Hut kuvvet hesabında, her bir parçacık ağaç düğümlerini gezer. Bu gezinim sırasında **bellek erişim kalıbı düzensizdir** (warp divergence). CUDA Tile şu avantajları sağlar:

1. **Yapılandırılmış Paylaşımlı Bellek:** Parçacık verileri tile'lara yüklenerek L1 cache miss oranı düşürülür
2. **Komşu Parçacık Gruplandırması:** Aynı octree yaprağındaki parçacıklar aynı tile'a yerleştirilerek veri yeniden kullanımı sağlanır
3. **Warp Divergence Azaltma:** Tile tabanlı erişim, koşullu dallanma yerine yapılandırılmış okuma kullanır

### Nasıl Uygulanır?

```
Mevcut durum (CPU):
  for each particle p:
    traverse tree → interact(p, node)

CUDA Tile yaklaşımı:
  1. Parçacıkları Morton kodu ile sırala (spatial locality)
  2. Her tile = 32 komşu parçacık
  3. Tile'lar shared memory'ye yüklenir
  4. Ağaç düğümleri ile etkileşim tile düzeyinde yapılır
```

### Uygulama Zorluğu: YÜKSEK
- Mevcut pointer-tabanlı ağaç yapısı GPU uyumlu değil
- Önce Structure-of-Arrays (SoA) dönüşümü gerekli
- Morton kodu sıralaması gerekli

### Beklenen Performans Etkisi
- Kuvvet hesabında **3-8x** hızlanma (düzensiz dağılımda daha az, düzenli dağılımda daha çok)
- Shared memory bant genişliği: Global memory'ye göre ~10x daha hızlı

---

## 3. Özellik 2: CCCL 3.2 Modern C++ API'leri

### Ne?
CCCL 3.2 (CUDA C++ Core Libraries), idiomatik C++ arayüzleri sunar:
- `cuda::stream` — RAII stream yönetimi
- `cuda::event` — Modern event senkronizasyonu
- `cuda::launch` — Tip-güvenli kernel başlatma
- Memory Resources — Özel bellek yönetimi

### Neden Bu Projeye Uygulanmalı?

Mevcut kod zaten C++20 kullanmaktadır. Ancak CUDA tarafı hâlâ C-style API kullanmaktadır:

```cpp
// Mevcut (C-style):
cudaStream_t stream;
cudaStreamCreate(&stream);
cudaMalloc(&d_particles_, count * sizeof(Particle));
cudaMemcpyAsync(d_particles_, particles, ...);
// ... hata kontrolü manuel

// CCCL 3.2 (Modern C++):
auto stream = cuda::stream{};
auto d_particles = cuda::memory::async_resource{stream}.allocate<Particle>(count);
cuda::launch(kernel, grid, block, stream, d_particles, count);
// RAII: otomatik temizlik, exception-safe
```

### Nasıl Uygulanır?

1. `cuda_renderer.cu` → CCCL 3.2 API'lerine taşıma
2. Yeni `bh_cuda_kernels.cu` → Çekirdek algoritma GPU kernel'leri
3. Benchmark framework → `cuda::stream` ve `cuda::event` kullanımı

### Uygulama Zorluğu: ORTA
- API değişikliği doğrudan (1:1 karşılık)
- Mevcut C++20 altyapısı ile uyumlu
- Hata yönetimi RAII ile basitleşir

### Beklenen Performans Etkisi
- Doğrudan performans etkisi düşük (API overhead azaltma)
- **Geliştirme verimliliği** açısından büyük kazanım
- Memory leak riskini ortadan kaldırır

---

## 4. Özellik 3: cub::DeviceTopK ve Segmented Reduction

### Ne?
CCCL 3.2 yeni algoritmalar sunar:
- `cub::DeviceTopK` — K-seçimi, radix sort'a göre **5x** hızlı
- `cub::DeviceSegmentedReduce` (sabit boyut) — Küçük segmentlerde **66x** hızlanma
- `cub::DeviceSegmentedScan` — Segment bazlı prefix scan
- `cub::DeviceFind` — İkili arama ve koşullu arama

### Neden Bu Projeye Uygulanmalı?

**1. Kütle Merkezi Hesabı (Upward Pass):**
Octree'nin her düğümü için kütle merkezi hesaplanır. Bu, segment bazlı reduction'dır:
- Her yaprak düğümü = 1 segment
- Her segment = düğümdeki parçacıkların kütle ve konum toplamı
- `cub::DeviceSegmentedReduce` bu işlemi **massively parallel** yapar

```
Mevcut (CPU, Recursive):
  compute_center_of_mass(root)
    → for each child: recursive call
    → sum masses, compute weighted average

CUDA 13.2 (cub::DeviceSegmentedReduce):
  1. Parçacıkları düğüm bazında segment'lere ayır
  2. Her segment için mass * position toplamını paralel hesapla
  3. Bottom-up level-by-level reduction
  → 66x hızlanma (küçük segment boyutları için)
```

**2. En Yakın Komşu / En Büyük Kuvvet Tespiti:**
`cub::DeviceTopK` ile en büyük kuvvet etkileşimini yaşayan K parçacığı hızlıca bulunabilir.

**3. Bounding Box Hesabı:**
`cub::DeviceReduce` ile min/max pozisyon hesabı tek kernel çağrısı ile yapılır.

### Uygulama Zorluğu: DÜŞÜK
- cub API'leri doğrudan uygulanabilir
- SoA veri yapısı gerektirir (position, mass ayrı diziler)
- Segment bilgisi (offset dizisi) hazırlanmalı

### Beklenen Performans Etkisi
- Upward pass: **10-66x** hızlanma (segment boyutuna bağlı)
- Bounding box: **5-20x** hızlanma
- TopK: Analiz/debug amaçlı, simülasyon hızını doğrudan etkilemez

---

## 5. Özellik 4: Geliştirilmiş Bellek Yönetimi API'leri

### Ne?
- `cudaMemcpyWithAttributesAsync` — Öznitelik tabanlı esnek bellek transferi
- `cudaMemPoolGetAttribute` — Bellek havuzu sorgulaması
- Per-context yerel bellek ayak izi azaltma

### Neden Bu Projeye Uygulanmalı?

Mevcut kodda her simülasyon adımında Host→Device transfer yapılmaktadır:

```cpp
// cuda_renderer.cu:266 — Her frame'de tam kopyalama
cudaMemcpyAsync(d_particles_, particles, count * sizeof(Particle),
                cudaMemcpyHostToDevice, stream_);
```

Yeni API'ler şunları sağlar:

1. **Kısmi Güncelleme:** Yalnızca değişen parçacık verilerini transfer et (pozisyon + hız, kuvvet hariç)
2. **Bellek Havuzu Yönetimi:** Simülasyon boyunca tek seferlik tahsis, sonra yeniden kullanım
3. **Öznitelik Tabanlı Transfer:** Önbellek ipuçları ile transfer optimizasyonu

### Nasıl Uygulanır?

```cpp
// Parçacık verilerinin yalnızca güncellenen kısımlarını transfer et
cudaMemcpyWithAttributesAsync(d_positions, h_positions,
    count * sizeof(float4), cudaMemcpyHostToDevice,
    {.srcAccessOrder = cudaMemAccessOrderRelaxed}, stream);
```

### Uygulama Zorluğu: DÜŞÜK
### Beklenen Performans Etkisi
- Transfer süresi: **%20-40** azalma
- Büyük parçacık sayılarında (>100K) belirgin fark

---

## 6. Özellik 5: C++20 Uyumluluk İyileştirmeleri (nvcc)

### Ne?
CUDA 13.2'de nvcc derleyicisi aşağıdaki C++20 düzeltmelerini içerir:
- Constraints ve requires-expressions düzeltmeleri
- Lambda ifadelerinde iyileştirmeler
- `noexcept` specification düzeltmeleri
- `[[no_unique_address]]` attribute düzeltmeleri
- Şablonlu ve miras alınmış kodda doğruluk iyileştirmeleri

### Neden Bu Projeye Uygulanmalı?

Mevcut proje tam C++20 kullanmaktadır ancak CUDA derlemesi C++17'de (`CMAKE_CUDA_STANDARD 17`). Bunun nedeni eski nvcc'nin C++20 desteğinin eksik olmasıydı.

CUDA 13.2 ile:
1. **`Vector3D` sınıfı doğrudan `__device__` olarak kullanılabilir** — `constexpr` operatörleri GPU'da çalışır
2. **`Particle` sınıfı GPU uyumlu** — `[[nodiscard]]`, `noexcept` doğru derlenir
3. **Concepts kullanımı** — GPU kernel'lerinde tip kısıtlamaları

```cpp
// Artık doğrudan GPU'da çalışır:
__device__ constexpr Vector3D operator+(const Vector3D& other) const noexcept {
    return Vector3D{data_[0]+other.data_[0], data_[1]+other.data_[1], data_[2]+other.data_[2]};
}
```

### Nasıl Uygulanır?

1. `CMAKE_CUDA_STANDARD` → 20 olarak güncelle
2. `vektor.h` ve `particle.h` → `__host__ __device__` qualifier ekle
3. GPU kernel'leri mevcut C++ sınıflarını doğrudan kullansın

### Uygulama Zorluğu: DÜŞÜK
### Beklenen Performans Etkisi
- Doğrudan performans etkisi yok
- **Kod tekrarını ortadan kaldırır** (GPU için ayrı veri yapıları gerekmez)
- Bakım kolaylığı büyük ölçüde artar

---

## 7. Özellik 6: CUDA Graphs Geliştirmeleri

### Ne?
CUDA 13.2 (cuda.core 0.6) ile CUDA Graphs:
- Koşullu yürütme (conditional execution)
- Fork-join kalıpları
- Polimorfik `cudaGraphNodeGetParams` API

### Neden Bu Projeye Uygulanmalı?

Barnes-Hut simülasyon döngüsü her adımda aynı kernel dizisini çalıştırır:

```
Her simülasyon adımı:
  1. find_bounding_box()      → GPU kernel
  2. build_tree()             → GPU kernel (veya CPU+transfer)
  3. compute_mass_distribution() → GPU kernel (segmented reduce)
  4. calculate_forces()       → GPU kernel (ana hesaplama)
  5. integrate_particles()    → GPU kernel
```

CUDA Graph ile bu 5 adımlık dizi **tek seferlik kaydedilir** ve her adımda yeniden başlatma overhead'i olmadan tekrar yürütülür.

### Nasıl Uygulanır?

```cpp
// Graph oluştur (bir kez)
cudaGraph_t graph;
cudaGraphExec_t graphExec;
cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);
  // Tüm kernel'leri başlat
  launch_bounding_box_kernel(...);
  launch_force_kernel(...);
  launch_integrate_kernel(...);
cudaStreamEndCapture(stream, &graph);
cudaGraphInstantiate(&graphExec, graph);

// Her adımda çalıştır (çok düşük overhead)
while (simulating) {
    cudaGraphLaunch(graphExec, stream);
    cudaStreamSynchronize(stream);
}
```

### Uygulama Zorluğu: ORTA
- Kernel bağımlılıkları doğru tanımlanmalı
- Dinamik ağaç yapısı Graph ile uyumsuz olabilir → sabit boyutlu yaklaşım gerekli

### Beklenen Performans Etkisi
- Kernel başlatma overhead'i: **%40-60** azalma
- Küçük parçacık sayılarında (<10K) çok belirgin
- Büyük sayılarda hesaplama dominant olduğu için daha az etkili

---

## 8. Özellik 7: Host Task Spin-Wait Modu

### Ne?
`cudaLaunchHostFunc()` ve graph host node'ları artık spin-wait dispatch modunu destekler. Bu, interrupt-based blocking'e göre düşük gecikme sağlar.

### Neden Bu Projeye Uygulanmalı?

Simülasyon döngüsünde CPU-GPU senkronizasyonu gereklidir:
- Ağaç CPU'da oluşturulur → GPU'ya transfer → GPU hesaplama → sonuç geri al

Spin-wait modu ile:
- Host callback'lerin gecikmesi **mikrosaniye** düzeyine iner
- CPU-GPU pipeline'ı daha sıkı bağlanır

### Uygulama Zorluğu: DÜŞÜK
### Beklenen Performans Etkisi
- Senkronizasyon gecikmesi: **%50-70** azalma
- Genel simülasyon etkisi: **%2-5** (senkronizasyon toplam sürenin küçük bir kısmı)

---

## 9. Özellik 8: Matematik Kütüphanesi Optimizasyonları

### Ne?
- `expm1f()`: %20'ye kadar hızlanma
- `erff()`: %5-10 hızlanma

### Neden Bu Projeye Uygulanmalı?

Barnes-Hut simülasyonunda temel matematik işlemleri `sqrt`, `pow`, çarpma ve bölme'dir. `expm1f` ve `erff` doğrudan kullanılmamaktadır.

**Dolaylı Etki:** Eğer simülasyona aşağıdakiler eklenirse:
- Termal dağılım (Maxwell-Boltzmann) → `erff` kullanır
- Enerji koruma analizi → `expm1f` kullanılabilir

### Uygulama Zorluğu: DÜŞÜK (ancak doğrudan fayda yok)
### Beklenen Performans Etkisi
- Mevcut simülasyon: **etki yok**
- Genişletilmiş fizik modeli ile: **%5-10** (ilgili hesaplamalarda)

---

## 10. Özellik 9: Gelişmiş Mimari Desteği (Blackwell)

### Ne?
CUDA 13.2, Compute Capability 10.x ve 12.x (Blackwell) mimarilerini tam destekler:
- cuBLAS'ta MXFP8 Grouped GEMM
- cuFFT'de power-of-2 optimizasyonları
- cuSOLVER'da FP64 fixed-point emulation

### Neden Bu Projeye Uygulanmalı?

1. **Blackwell GPU'larda** kuvvet hesabı kernel'leri daha verimli çalışır
2. **Tensor Core kullanımı** (MXFP8) → N-body simülasyonunda matris çarpımı ile hızlı kuvvet hesabı
3. CMakeLists.txt'e `100` (Blackwell) mimarisi eklenmeli

### Nasıl Uygulanır?

```cmake
# CMakeLists.txt güncellemesi
set(CMAKE_CUDA_ARCHITECTURES 75 80 86 89 100 120)  # + Blackwell
```

### Uygulama Zorluğu: DÜŞÜK
### Beklenen Performans Etkisi
- Blackwell donanımı varsa: **%20-50** genel hızlanma (yeni SM mimarisi)
- Mevcut donanımda: etki yok

---

## 11. Özellik 10: Birleşik ARM Toolkit

### Ne?
CUDA 13.2'den itibaren aynı Arm SBSA CUDA Toolkit, sunucu ve gömülü (Jetson Thor/Orin) cihazlarda kullanılabilir.

### Neden Bu Projeye Uygulanmalı?

- Jetson platformlarında N-body simülasyonu çalıştırma imkanı
- CI/CD pipeline basitleştirme (tek toolkit, çoklu hedef)
- Edge computing uygulamaları (ör. uydu yörünge simülasyonu gömülü sistemde)

### Uygulama Zorluğu: DÜŞÜK
### Beklenen Performans Etkisi
- ARM platformlarında: Tam GPU hızlandırma desteği
- x86 platformlarında: Etki yok

---

## 12. Geçiş Stratejisi ve Öncelik Matrisi

### Faz 1: Temel Altyapı (1-2 hafta)
1. CMakeLists.txt → CUDA 13.2, C++20 CUDA standardı
2. Veri yapılarını SoA (Structure-of-Arrays) formatına dönüştür
3. CCCL 3.2 Modern API'lerine geçiş

### Faz 2: Çekirdek GPU Kernel'leri (2-3 hafta)
4. Bounding box → `cub::DeviceReduce` ile GPU'ya taşı
5. Kütle merkezi → `cub::DeviceSegmentedReduce` ile GPU'ya taşı
6. Kuvvet hesabı → CUDA kernel (tile-optimized)
7. Integrasyon → Basit CUDA kernel

### Faz 3: Optimizasyon (1-2 hafta)
8. CUDA Graphs ile simülasyon döngüsü
9. Bellek yönetimi optimizasyonu
10. Blackwell mimari desteği

### Faz 4: Test ve Doğrulama (1 hafta)
11. CPU vs GPU doğruluk karşılaştırması
12. Performans benchmark'ları
13. Düzenli/düzensiz dağılım testleri

---

## 13. Karşılaştırmalı Performans Tahmin Tablosu

| Bileşen | CPU Seri | CPU OpenMP (8T) | GPU (Mevcut) | GPU (CUDA 13.2) |
|---------|----------|-----------------|---------------|------------------|
| Bounding Box | 1x | ~4x | N/A | ~20x (cub Reduce) |
| Ağaç Oluşturma | 1x | 1x (seri) | N/A | ~3x (sort+build) |
| Kütle Merkezi | 1x | ~4x | N/A | ~30x (seg. reduce) |
| Kuvvet Hesabı | 1x | ~6x | N/A | ~15-40x (tile kernel) |
| İntegrasyon | 1x | ~7x | N/A | ~50x (trivial parallel) |
| **Toplam Adım** | **1x** | **~5x** | **N/A** | **~12-30x** |

> **Not:** GPU performansı parçacık dağılımına güçlü şekilde bağımlıdır.
> Düzensiz dağılımda warp divergence nedeniyle GPU avantajı azalır.

### Parçacık Sayısına Göre Beklenen Hızlanma

| N (Parçacık) | CPU Seri | CPU OpenMP | GPU CUDA 13.2 |
|-------------|----------|------------|---------------|
| 1,000 | 1x | 3x | 2x (overhead dominant) |
| 10,000 | 1x | 5x | 8x |
| 100,000 | 1x | 6x | 25x |
| 1,000,000 | 1x | 6x | 40x |
| 10,000,000 | 1x | 6.5x | 60x+ |

---

## Kaynaklar

- [CUDA 13.2 Release Notes](https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html)
- [CUDA 13.2 Blog Post](https://developer.nvidia.com/blog/cuda-13-2-introduces-enhanced-cuda-tile-support-and-new-python-features/)
- [CCCL 3.2 Documentation](https://nvidia.github.io/cccl/)
- [Barnes-Hut GPU Algorithms — Burtscher & Pingali (2011)](https://iss.oden.utexas.edu/Publications/Papers/burtscher11.pdf)
