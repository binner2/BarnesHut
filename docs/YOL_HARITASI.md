# Yol Haritası ve Durum

Güncelleme: 5 Kasım 2026. Dal: `claude/sweet-faraday-fyldry`.

## Ne yapıldı, ne yapılmadı

| Soru | Durum | Kanıt |
|------|-------|-------|
| Kodlar düzeltildi mi? | **Evet.** Depo derlenmiyordu; 4 hata giderildi (`generate_random` bildirimi, `omp.h`, `quadtree_visualizer.h` varsayılan argüman ve `length()`, CMake'te zorunlu CUDA). CUDA/OpenGL isteğe bağlı. | `cmake --build` temiz, `ctest` 1/1 geçiyor |
| Nokta bulutu çözümü yazıldı mı? | **Evet.** `pointcloud/` modülü: IO, lineer octree, BH-KDE, BH-mean-shift, BH-WLOP, normal, SOR, voxel, CLI, 6 test grubu | `pointcloud/README.md` |
| Simülasyon deneyleri hazırlandı mı? | **Kısmen → şimdi evet.** Sentetik sahne, benchmark komutu ve 4 Kasım ölçümleri vardı; tekrar üretilebilir betikler bu güncellemeyle eklendi | `experiments/run_benchmarks.sh`, `experiments/plot_benchmarks.py`, `experiments/results/2026-11-04/` |
| Makale yazıldı mı? | **Hayır.** Yazılan şey teknik rapor + araştırma planıdır (`docs/derinlik_kamerasi_barnes_hut_coworks.html`). Makale taslağı için gerçek veri sonuçları ve taban çizgileri gerekir | — |
| Gerçek veri planı ayrı dosyada mı? | **Yoktu → şimdi var.** | `docs/GERCEK_VERI_PLANI.md` |

## Sırada yapılacaklar (öncelik sırasıyla)

1. **Gerçek veri** — `docs/GERCEK_VERI_PLANI.md` bölüm 0–1: TUM `freiburg1_xyz` ile başla; aynı CSV şemasıyla kare başına süre ve hata. (1–2 hafta)
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
```

Rapor kaynağı `docs/report_src/` altındadır; `assemble_report.py` CSV'lerden
sayıları, tabloları ve SVG grafiklerini şablona işler, ardından `premium-html`
build betiğiyle HTML üretilir.
