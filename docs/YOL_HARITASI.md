# Yol Haritası ve Durum

Güncelleme: 13 Kasım 2026. Dal: `ccr-aef4f1eb-5e6s61` (`claude/sweet-faraday-fyldry` üzerine).

## Ne yapıldı, ne yapılmadı

| Soru | Durum | Kanıt |
|------|-------|-------|
| Kodlar düzeltildi mi? | **Evet.** Depo derlenmiyordu; 4 hata giderildi (`generate_random` bildirimi, `omp.h`, `quadtree_visualizer.h` varsayılan argüman ve `length()`, CMake'te zorunlu CUDA). CUDA/OpenGL isteğe bağlı. | `cmake --build` temiz, `ctest` 1/1 geçiyor |
| Nokta bulutu çözümü yazıldı mı? | **Evet.** `pointcloud/` modülü: IO, lineer octree, BH-KDE, BH-mean-shift, BH-WLOP, normal, SOR, voxel, CLI, 6 test grubu | `pointcloud/README.md` |
| Simülasyon deneyleri hazırlandı mı? | **Kısmen → şimdi evet.** Sentetik sahne, benchmark komutu ve 4 Kasım ölçümleri vardı; tekrar üretilebilir betikler bu güncellemeyle eklendi | `experiments/run_benchmarks.sh`, `experiments/plot_benchmarks.py`, `experiments/results/2026-11-04/` |
| Makale yazıldı mı? | **Hayır.** Yazılan şey teknik rapor + araştırma planıdır (`docs/derinlik_kamerasi_barnes_hut_coworks.html`, sürüm 1.1). Makale taslağı için taban çizgileri ve ek diziler gerekir | — |
| Gerçek veri planı ayrı dosyada mı? | **Evet.** Bölüm 0, 1.1 ve 3 uygulandı; 1.2–1.4, 2 ve 4'ün çoğu bekliyor | `docs/GERCEK_VERI_PLANI.md` |
| Gerçek veri ölçüldü mü? | **Kısmen → evet (tek dizi).** TUM `freiburg1_xyz`, 20 eşit aralıklı kare, stride 2, 5 tekrar medyanı: octree inşası, KDE kesin/θ/κ (h = 2 ve 5 cm), mean shift θ∈{0, 0,5, 0,7, 1,0}, WLOP θ∈{0, 0,5, 0,7, 1,0}; tam çözünürlükte derinlik–NN profili. Sentetik üç bulgu (κ kazancı, θ = 1,0 bozulması, 1/z² yoğunluk) gerçek veride doğrulandı; WLOP gerçek karede CoV'yi artırdı (neden ölçülmedi) | `experiments/results/2026-11-13/`, rapor "Gerçek veri" bölümü |

## Sırada yapılacaklar (öncelik sırasıyla)

1. **Gerçek veri, devam** — `freiburg1_xyz` bitti. Sırada: `freiburg1_desk`, `freiburg2_xyz`, `freiburg3_*` (aynı `experiments/run_tum.sh`); kalibre iç parametre farkı; WLOP'un gerçek karede CoV'yi artırmasının nedeni (yoğunluk ağırlığı, yineleme sayısı); mean shift erken durdurma. (1 hafta)
2. **Taban çizgileri** — Open3D/PCL yarıçap KDE, Poisson-disk, CGAL WLOP, IFGT; aynı parametrelerle. (1 hafta)
3. **Eksen A'nın teorisi** — ikinci mertebe kestirimi kesin üst sınıra çevir (üçüncü mertebe kalan), \(\theta\)–\(\kappa\)–yaprak için seçim kuralı; sentetik + TUM Pareto eğrileri. (2 hafta)
4. **Kod eksikleri** — `docs/GERCEK_VERI_PLANI.md` bölüm 4 (değişken \(h(z)\), 6-D mean shift, erken durdurma, PCD compressed). (2 hafta)
5. **Makale taslağı** — IEEE konferans biçimi, eksen A (+B): giriş, ilgili çalışma (sistematik yakın-rakip taraması yapılacak), yöntem (rapor bölüm 2–4'ten), deneyler (sentetik + gerçek), tartışma. Taslağa başlamadan önce 1–3 bitmiş olmalı; aksi hâlde deney bölümü boş kalır.
6. **GPU (eksen D)** — Karras tipi inşa + warp-uyumlu `traverse_bh`; mevcut `visualization/cuda_renderer.cu` ile birleştirme. (4+ hafta, CUDA makinesi gerekir)

## Tekrar üretim

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j && ctest --test-dir build
experiments/run_benchmarks.sh build/pointcloud/pointcloud_tool experiments/results/$(date +%F)
python3 experiments/plot_benchmarks.py experiments/results/$(date +%F) experiments/results/$(date +%F)/fig
# gerçek veri (TUM dizisi data/tum/ altına açılmış olmalı; ~2 saat, 4 iş parçacığı)
experiments/run_tum.sh build/pointcloud/pointcloud_tool data/tum/rgbd_dataset_freiburg1_xyz experiments/results/$(date +%F) 20 2 5
python3 experiments/plot_tum.py experiments/results/$(date +%F) experiments/results/$(date +%F)/fig_tum
```

Rapor kaynağı `docs/report_src/` altındadır; `assemble_report.py` sentetik ve TUM
CSV'lerinden sayıları, tabloları ve SVG grafiklerini şablona işler, ardından
`premium-html` build betiğiyle HTML üretilir (`docs/report_src/README.md`).
