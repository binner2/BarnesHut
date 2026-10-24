# CUDA 13.4 modernizasyonu: denetim, değişiklikler, kabul ölçütleri

`cuda/` (eski `cuda_13_2/`) · CUDA 13.4 Update 1 (nvcc 13.4.92, CCCL 3.4.3) · 24 Ekim 2026

Bu belge, `claude/cuda-13-2-integration-dgXHI` branch'indeki kodu (bundan sonra **M0**, commit `0fbd396`) satır satır inceleyerek hazırlandı. Kodun gerçek durumunu, CUDA 13.4'ten gerçekte neyin alındığını, yapılan değişiklikleri ve her değişikliğin kabul ölçütünü listeler.

**Doğrulama kapsamı.** Kodun tamamı nvcc 13.4.92 ile derlendi: sm_75…sm_121 için SASS ve compute_121 için PTX üretildi. GPU algoritması `__host__ __device__` fonksiyonlardan oluştuğu için host'ta adım adım emüle edilerek test edildi (`cuda/tests/host_core_test.cu`, 72/72 test). CPU tarafı ASan/UBSan ve OpenMP altında test edildi. Ancak **bu ortamda GPU yok**. GPU üzerinde çalışma, zamanlama, CUDA Graph eşdeğerliği ve compute-sanitizer çalıştırması henüz doğrulanmadı ve bu belgede **[GPU'da doğrulanmadı]** olarak işaretlendi.

---

## 1. M0'ın etiketlediği "CUDA 13.2 özellikleri" gerçekte ne?

| M0'daki iddia (başlık yorumları, rapor) | Koddaki gerçek durum | Karar |
|---|---|---|
| "cub::DeviceReduce for bounding box" | Yalnızca yorumda geçiyor. Bbox, el yazımı shared-memory + `atomicMin/Max` kernel'i ile hesaplanıyor. Sonuç *ordered-int* biçiminde bırakılıyor ve hiçbir yerde float'a geri çevrilmiyor (`ordered_int_to_float` tanımlı ama çağrılmıyor). | Gerçekten `cub::DeviceReduce` + `cuda::transform_iterator` ile yeniden yazıldı |
| "cub::DeviceSegmentedReduce for mass distribution" | Yalnızca yorumda geçiyor. Kütle dağılımı, her ağaç seviyesi için tüm node'lar üzerinde ayrı bir kernel başlatıyor. | Kaldırıldı. Yerine tek geçişli bottom-up refit geldi |
| "Tile-optimized force calculation" | `force_kernel_tiled` tiling ya da shared memory kullanmıyor. **Hiçbir yerden çağrılmıyor.** | Kaldırıldı. Yerine yeni traversal kernel'i geldi |
| "CUDA Graphs for simulation loop" | M0'da `cudaGraph*` çağrısı yok. | Gerçek graph capture/replay eklendi (opsiyonel) |
| "Modern CCCL 3.2 API usage" | Tek CUB çağrısı olan `DeviceRadixSort::SortPairs` hatalı kullanılıyor: in-place (aynı giriş/çıkış buffer'ı), geçici bellek boyutu sorgulanmıyor ("conservative estimate"), dönüş kodu kontrol edilmiyor. Fonksiyon **hiç çağrılmıyor**. Sıralanan permütasyon da veriye uygulanmıyor. | Doğru two-phase kullanım + gather |
| `benchmark_gpu_bh` ("GPU Barnes-Hut") | `benchmark_gpu_direct`'i çağırıyor. GPU_BH olarak raporlanan her sayı aslında O(N²) direct-sum'a ait. | Gerçek LBVH Barnes–Hut yolu |

## 2. Bulgular (önem sırasıyla)

K: kritik (sonuç veya iddia yanlış) · O: optimizasyon · M: ölçüm/derleme

| # | Yer (M0) | Bulgu | Sonuç | Durum |
|---|---|---|---|---|
| K1 | `benchmark_gpu_bh` | GPU BH, direct-sum'a yönlendiriliyor | GPU BH için raporlanmış hiçbir sayı geçerli değil | Düzeltildi |
| K2 | CPU çekirdeği (`tree.cpp`, `particle.h`) | Node'lar hem pool hem `unique_ptr` children tarafından sahipleniliyor | İkinci adımda double free. Benchmark'ın CPU referansı çöküyor (ASan: heap-use-after-free) | Düzeltildi |
| K3 | `CORE_SOURCES` | `file.cpp` (`generate_random` bildirimi yok) ve `BHtreetest.cpp` (`omp.h` yok) derlenmiyor | `benchmark_cuda13` hedefi M0'da **hiç link edilemiyordu** | Düzeltildi |
| K4 | `integrate_kernel`, `Particle::integrate` | "KDK leapfrog" iki yarım kick'te de aynı \(a(x_n)\)'yi kullanıyor: hız için birinci mertebe Euler | İki-cisim yörüngesinde 20 000 adımda \(\lvert\Delta E/E\rvert\) = 0.58 | Düzeltildi: 1e-8 (CPU, double), 2e-5 (GPU çekirdeği, float) |
| K5 | `force_kernel_tiled` | Stack taştığında node'lar **sessizce atlanıyor** | Hatalı kuvvet, hata bildirimi yok | Yeni kernel'de taşma monopole ile ele alınıyor ve sayılıyor (`stack_overflows`) |
| K6 | Kabul kriteri (CPU ve GPU) | Parçacığı içeren hücre kabul edilebiliyor | θ ≳ 0.58'de self-force | Düzeltildi (inside-box testi) |
| K7 | `benchmark_cpu_serial` | Hızlar geri kopyalanmıyor | Enerji testi sabit kinetik enerjiyle ölçüyordu | Düzeltildi |
| K8 | `run_benchmark_suite` | Doğruluk, K adım sonra *farklı* yörüngelerdeki kuvvetlerle ve BH (yaklaşık) referansla ölçülüyor | Doğruluk sütunu anlamsız | Ortak başlangıç durumunda double direct-sum'a karşı ölçülüyor |
| K9 | `--correctness-only` | `test_gpu_correctness` bu modda hiç çağrılmıyor | GPU doğruluğu CI'da test edilmiyordu | Düzeltildi |
| K10 | CPU OpenMP | `stats_` sayaçlarına tüm thread'ler yazıyor | Data race + cache-line contention: 4 thread'de 0.98× | `reduction`: 4.03× (N = 200k) |
| K11 | Çakışık parçacıklar | Maksimum ağaç derinliği yok | `test_edge_cases` (aynı konumda 100 parçacık) ağaç inşasında sonlanmıyor (ownership düzeltilmiş PR #6 kodunda 60 s timeout ile doğrulandı) | `MAX_TREE_DEPTH` |
| M1 | Benchmark | Warmup adımları da zamanlanıyor, mean ± std raporlanıyor | Isınma etkisi ölçüme karışıyor | Warmup zamanlanmıyor; adım başına örnek, medyan + IQR |
| M2 | `CPU_Serial` | OpenMP thread'leriyle çalışıyor (`CPU_OpenMP` ile aynı fonksiyon) | "OpenMP speedup" ≈ 1× (tanım gereği) | Serial artık 1 thread |
| M3 | CMake | `CUDA_SEPARABLE_COMPILATION` (rdc) gereksiz; `-Wall`/`-march=native` nvcc'ye de geçiyor; `-ffast-math` NaN/Inf testlerini geçersiz kılıyor; sm_121 (GB10) listede yok; OpenGL `REQUIRED` | — | Düzeltildi (bkz. §5) |
| M4 | Zamanlama | Her adımda 5 event oluşturulup yok ediliyor | Gereksiz host overhead'i | Event'ler bir kez oluşturuluyor |
| O1 | Mass distribution | Seviye sayısı kadar kernel, her biri tüm node'ları tarıyor | O(N·depth) iş | Tek kernel, O(N) |
| O2 | Bbox / Morton | Bbox host'a okunmadan Morton kernel'ine geçiyor (M0'da bbox hiç okunmuyordu) | — | Bbox device'ta kalıyor, kernel doğrudan okuyor (graph-uyumlu) |

## 3. CUDA 13.4'ten projeye gerçekten alınanlar

Kaynak: CUDA 13.4 / 13.4 Update 1 release notes. Sürüm numaraları kurulu toolkit ile doğrulandı: `nvcc --version` → V13.4.92, `CCCL_VERSION` → 3004003.

| 13.4 öğesi | Projede nereye | Karar |
|---|---|---|
| CCCL 3.4: Thrust iterator'ları `cuda::` eşdeğerleri lehine deprecated | Bbox girişinde `cuda::counting_iterator` + `cuda::make_transform_iterator`. Ara dizi oluşmuyor. | **UYGULANDI** |
| CCCL 3.4: tüm device algoritmaları için tek çağrılı, environment tabanlı overload'lar (`EnvT env = {}`) | Header'da mevcut (`device_radix_sort.cuh`) | **UYGULANMADI.** Graph capture için geçici bellek önceden ve sabit boyutta ayrılmalı, bu yüzden two-phase API bilinçli olarak korundu |
| CCCL 3.4: `cuda::std` parallel algoritmaları (`cuda::execution::gpu`) | Gather/bbox için aday | Ertelendi (runtime davranışı GPU'da ölçülmeden alınmadı) |
| CCCL 3.4: `cub::WarpReduceBatched` | Uygun bir iş yükü yok | UYGULANMADI |
| `-fdevice-sanitize=memcheck` (nvcc 13.4 `--help` ile doğrulandı) | `-DBH_CUDA_SANITIZE=ON` build seçeneği | **UYGULANDI** (derleme doğrulandı; çalıştırma [GPU'da doğrulanmadı]) |
| Nsight Compute 2026.3 | Ölçüm protokolü (§6) | ARAÇ |
| PTX ISA 9.4 | Kernel'ler ISA'ya özgü komut kullanmıyor | Yalnızca not |
| Rubin `sm_107` (developer preview) | `nvcc --list-gpu-arch` içinde var | Varsayılan fatbin'e EKLENMEDİ (preview) |
| CUDA Tile / cuTile Python 1.5–1.6, Tile IR | İş yükü düzensiz ağaç traversal'ı; yoğun tile işlemi yok | UYGULANMADI |
| Green contexts'in graph node'larıyla ilişkilendirilmesi | Ölçümde SM sayısını değiştirir | UYGULANMADI |
| MPS v3, Compute Fabric Transport, `cudaMemGetLocationInfo` | Tek GPU, tek süreç | Kapsam dışı |
| Sürücü: yeni özellikler için R615; 13.x minor-version uyumluluğu için ≥ 580; Linux'ta sürücü artık toolkit ile gelmiyor | Kod 13.4'e özgü bir runtime API'si kullanmıyor, ≥ 580 yeterli olmalı | [GPU'da doğrulanmadı] |

Kernel'lerin kullandığı hiçbir özellik 13.4'e özgü değildir: `cuda::atomic_ref` (libcu++), CUB, stream capture ve CUDA Graphs uzun süredir mevcuttur. Bu nedenle kod ve makale metninde "CUDA 13.4 özellikleri" ifadesi kullanılmamalıdır. Doğru ifade şudur: "CUDA 13.4 araç zinciri (nvcc 13.4.92, CCCL 3.4) ile derlenmiş, mimariye özgü özellik kullanmayan SIMT kodu".

## 4. Yeni GPU hattı

Her kuvvet değerlendirmesi tamamen device üzerinde çalışır. Host'a veri dönüşü olmadığı için bütün adım graph'a capture edilebilir:

1. **Bbox:** `cub::DeviceReduce::Reduce(transform_iterator(counting_iterator, PointBox), BoxMerge)`
2. **Morton anahtarları:** 63 bit (21 bit/eksen, magic-bits `split3_21`), kübik domain
3. **Sıralama:** `cub::DeviceRadixSort::SortPairs`, bitler [0, 63), ayrı giriş/çıkış buffer'ları
4. **Gather:** parçacık durumu Morton sırasına permüte edilir ve adımlar arasında bu sırada kalır. `download()` orijinal sırayı geri kurar.
5. **Hiyerarşi:** Karras (2012) binary radix tree. Her internal node için bir thread. Eşit anahtarlar indeksle ayrıştırılır (çakışık parçacıklar).
6. **Refit:** yapraktan köke doğru, `cuda::atomic_ref<int, thread_scope_device>::fetch_add(acq_rel)`. İlk gelen thread durur, ikinci gelen birleştirir. PTX'te `atom.add.acq_rel.gpu.s32` üretildi; node verisi coherent `ld.global` ile okunuyor (`.nc` değil).
7. **Kuvvet:** parçacık başına stack'li DFS (64 giriş). Kabul koşulu \(d^2 > (l/\theta)^2\) ve parçacığın node kutusunun dışında olması (\(l\): sıkı kutunun en uzun kenarı). Taşma olursa monopole kullanılır ve sayılır.
8. **Zaman integrasyonu:** KDK, `kick_drift_kernel` (birleşik) → yeniden ağaç + kuvvet → `kick_kernel`.

### Kabul kriteri seçimi (ölçüm)

Host emülasyonunda (N = 4000, ε² = 1e-4), aynı maliyette (parçacık başına P–P + P–C etkileşim) medyan bağıl kuvvet hatası karşılaştırıldı. `l` = sıkı kutunun en uzun kenarı; \(b_{max}\) = kütle merkezinden en uzak köşeye mesafe.

| Dağılım | `l`-kriteri: maliyet → medyan | \(b_{max}\)-kriteri: maliyet → medyan |
|---|---|---|
| Uniform | 314.8 → 4.92e-3 | 341.9 → 5.54e-3 |
| Plummer | 460.1 → 4.80e-3 | 468.3 → 7.40e-3 |
| Clustered | 233.8 → 5.52e-3 | 245.7 → 6.49e-3 |

`l`-kriteri üç dağılımda da aynı maliyette daha doğru. \(b_{max}\)'ın sağladığı teorik hata sınırı bu dağılımlarda pratik bir kazanca dönüşmedi. Self-force riski inside-box testiyle ayrıca kapatıldı.

### Ağaç yapısının doğrulanması (host emülasyonu)

4 dağılım × N ∈ {2, 3, 17, 1000, 50 000} için şunlar kontrol edildi: her node'a tam bir kez ulaşılıyor, parent/child tutarlı, kök kütlesi ve kütle merkezi double referansla eşleşiyor, her kutu çocuklarının kutularını kapsıyor. Ölçülen en büyük derinlik 33 (Plummer, N = 50 000); stack = 64.

### Kuvvet doğruluğu (host emülasyonu, N = 4000, float)

| Dağılım | θ = 0 (medyan / max) | θ = 0.5 (medyan / p99) | θ = 0.7 (medyan / p99) |
|---|---|---|---|
| Uniform | 9.1e-7 / 1.1e-5 | 8.4e-3 / 5.1e-2 | 1.9e-2 / 1.1e-1 |
| Plummer | 8.4e-7 / 4.7e-5 | 4.8e-3 / 4.3e-2 | 1.0e-2 / 9.9e-2 |
| Clustered | 3.3e-6 / 3.5e-5 | 9.5e-3 / 6.7e-2 | 2.1e-2 / 1.5e-1 |

Not: θ'nın anlamı CPU octree'den farklıdır (küp kenarı yerine sıkı kutunun en uzun kenarı; binary ağaç). CPU ve GPU sonuçlarını aynı θ'da değil, aynı doğrulukta karşılaştırın.

## 5. Derleme

```bash
export PATH=/usr/local/cuda-13.4/bin:$PATH        # nvcc 13.4.92
cmake -S . -B build                                # varsayılan: 75..121 SASS + compute_121 PTX
cmake --build build -j
ctest --test-dir build --output-on-failure         # host_core (GPU gerekmez) + correctness

# Tek mimari (ör. RTX 4090 + GB10):  -DCMAKE_CUDA_ARCHITECTURES="89;121"
# Sanitizer build:                    -DBH_CUDA_SANITIZE=ON -DCMAKE_BUILD_TYPE=Debug
#   compute-sanitizer --tool memcheck build/benchmark_cuda --correctness-only
# Benchmark (kanonik protokol):
#   build/benchmark_cuda --benchmark --warmup 50 --steps 100 --n 10000,100000,1000000
```

Bu ortamda yapılan derleme denetimi (nvcc 13.4.92, `-O3 -lineinfo`):
- `cuobjdump --list-elf`: sm_75, 80, 86, 89, 90, 100, 120 ve 121 için birer cubin ile compute_121 PTX.
- `force_kernel`: 31–35 register, 256 B stack (traversal stack'i), 0 spill. Diğer tüm kernel'ler 0 spill.
- `-fdevice-sanitize=memcheck` ile derleme hatasız.

CMake değişiklikleri:
- CUDA opsiyonel; toolkit yoksa CPU hedefleri yine derleniyor.
- Minimum `CUDAToolkit 13.0`; 13.4'ten eski bir sürüm bulunursa uyarı veriliyor.
- `-rdc` kaldırıldı.
- Host bayrakları yalnızca CXX'e uygulanıyor.
- `-ffast-math` ve `-use_fast_math` artık opsiyonel (`BH_FAST_MATH`, varsayılan kapalı).
- OpenGL opsiyonel.
- `ctest` eklendi.

## 6. Ölçüm protokolü

- **Warmup:** ölçülmeyen adımlar (`--warmup`). **Ölçüm:** her adım ayrı bir örnek (`--steps`). Raporlanan değer medyan ve IQR (Q3 − Q1). Kanonik protokol: `--warmup 50 --steps 100`.
- **Zamanlama:** GPU'da CUDA event'leri; eager modda faz kırılımıyla (tree / force / integrate), graph modunda yalnızca toplam süre. CPU'da `steady_clock`.
- **Doğruluk:** ortak başlangıç durumundaki kuvvetler, double direct-sum'a karşı. Geçme ölçütü global bağıl L2 hatası \(\lVert F-F_{ref}\rVert/\lVert F_{ref}\rVert\). Parçacık başına maksimum hata raporlanır ama ölçüt olarak kullanılmaz, çünkü net kuvveti sıfıra yakın parçacıklar bu değeri domine eder.
- **Kapsam:** `GPU_Direct` yalnızca kuvvet değerlendirmesini ölçer. `GPU_BH` ve `GPU_BH_Graph` tam KDK adımını ölçer. Bu ikisi birbirine bölünerek karşılaştırılmamalıdır.
- **ncu:** `--clock-control base --cache-control all`. Traversal'daki stall nedenlerini (bağımlı yükler, divergence) PM sampling ile ayırın.

## 7. Fazlar ve kabul ölçütleri

| Varyant | İçerik | Kabul ölçütü | Durum |
|---|---|---|---|
| **M0** | `0fbd396` (orijinal CUDA branch) | — (referans) | Benchmark link edilemiyor; GPU BH = direct |
| **M1** | CPU düzeltmeleri, CMake, benchmark protokolü | ASan temiz; θ = 0'da CPU = direct (1e-15); OpenMP ölçekleniyor; enerji sınırlı | **Sağlandı** (bu ortamda) |
| **M2** | GPU LBVH Barnes–Hut (eager) | Host emülasyonunda 72/72; GPU'da `--correctness-only` geçer; GPU direct L2 < 1e-4; GPU BH L2 < 2e-2 (θ = 0.5) | Host: **sağlandı**. GPU: [GPU'da doğrulanmadı] |
| **M3** | CUDA Graph replay | 10 adım sonunda eager ile bit-bit aynı konum/hız (test mevcut); küçük N'de toplam süre azalır | [GPU'da doğrulanmadı] |
| M4 (öneri) | Warp-coherent traversal veya yaprak bucket'ları (≥ 8 parçacık/yaprak), stack'in shared memory'ye taşınması | M2 ile aynı doğruluk; ncu'da force kernel süresi düşer | Yapılmadı |
| M5 (öneri) | Quadrupole momentleri | Aynı doğruluk için daha büyük θ, daha az etkileşim | Yapılmadı |

Hipotez sütunu ölçüm değildir. Sayıları GPU'da M1–M3 için yeniden üretmeniz gerekir. M0'dan alınmış "GPU BH" sayıları (K1) hiçbir yerde kullanılmamalıdır.

## Kaynaklar

- CUDA Toolkit 13.4 Release Notes: https://docs.nvidia.com/cuda/archive/13.4.1/cuda-toolkit-release-notes/index.html
- CUDA Toolkit 13.4 Update 1 Release Notes: https://docs.nvidia.com/cuda/archive/13.4.2/cuda-toolkit-release-notes/index.html
- T. Karras, "Maximizing parallelism in the construction of BVHs, octrees, and k-d trees", HPG 2012. [DOĞRULANMADI: künye/DOI]
- J. Barnes, P. Hut, "A hierarchical O(N log N) force-calculation algorithm", Nature 324 (1986). [DOĞRULANMADI: DOI]
- J. Barnes, "A modified tree code: don't laugh; it runs", J. Comput. Phys. (1990) ve J. K. Salmon, M. S. Warren, "Skeletons from the treecode closet", J. Comput. Phys. 111 (1994), \(b_{max}\) kriteri için. [DOĞRULANMADI]
