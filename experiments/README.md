# Deneyler

- `run_benchmarks.sh <pointcloud_tool> <çıktı_dizini>`: raporun bütün sentetik
  ölçümlerini (N × θ × κ taraması, yaprak kapasitesi, uçtan uca demo) yeniden üretir
  ve makine bilgisini `machine.txt`'ye yazar.
- `plot_benchmarks.py <sonuç_dizini> <önek>`: CSV'den Pareto ve ölçeklenme SVG'leri.
- `render_views.py <girdi.ply> <çıktı.webp> [splat]`: renkli PLY'yi sentetik kamera
  iç parametreleriyle görüntüye yeniden izdüşürür (rapordaki 4 panel).
- `results/2026-11-04/`: raporda kullanılan ham sentetik ölçümler (4 iş parçacığı, GCC 13.3,
  -O3 -march=native, bulut konteyneri; mutlak süreler donanıma bağlıdır).
- `run_tum.sh <pointcloud_tool> <dizi_dizini> <çıktı_dizini> [kare=20] [stride=2] [tekrar=5]`:
  gerçek veri protokolü (`docs/GERCEK_VERI_PLANI.md` bölüm 3). TUM RGB-D dizisinin bütün
  derinlik PNG'lerini `tools/tum_depth_to_ply.sh` ile PLY'ye çevirir, eşit aralıklı kareleri
  seçer ve her kare için `pointcloud_tool info` (tam çözünürlük istatistikleri, derinlik–NN
  profili) ile `pointcloud_tool benchframe` (5 tekrarın medyanı) çalıştırır. Yarıda kesilirse
  aynı komutla kaldığı kareden devam eder.
- `plot_tum.py <sonuç_dizini> <önek> [--inline]`: TUM CSV'lerinden kare başına KDE hatası,
  mean shift etiket uyumu ve derinlik–NN profili SVG'leri.
- `results/2026-11-13/`: TUM `freiburg1_xyz` ölçümleri (`tum_freiburg1_xyz_frames.csv`,
  `tum_freiburg1_xyz_info.csv`, `tum_freiburg1_xyz_info_stride2.csv`, `profiles/`,
  `machine.txt`, günlükler). Benchmark bulutları stride 2 (her iki eksende her ikinci piksel,
  ~57 k nokta/kare); tam çözünürlük (~230 k) yalnız istatistik ve derinlik profili için.

CSV şeması: `algorithm,n,theta,kappa,cutoff_sigmas,leaf,threads,build_s,time_s,
far_evals,near_evals,max_rel_err,mean_rel_err,rms_rel_err,extra`.
`extra` alanı `anahtar=değer;...` biçimindedir (mean shift: `clusters`, `agreement`;
WLOP: `nn_cov`, `m`). `benchframe` satırlarında ek olarak `frame` (zaman damgası), `h`,
`repeats`, `p05`/`p95` (tekrarların 5–95 persentil süresi; `time_s` medyandır), mean shift
için `mean_iters`, WLOP için `nn_cov_input` ve `nn_mean`, `octree` satırında `nodes`, `depth`
bulunur. Bu satırlarda `n` ölçülen bulutun nokta sayısıdır.

`tum_*_info*.csv` şeması: `tag,n,min_x,min_y,min_z,max_x,max_y,max_z,nodes,leaves,depth,
root_width,build_s,nn_mean,nn_std,nn_cov,nn_min,nn_max`. `profiles/<kare>.csv` şeması:
`tag,z_lo,z_hi,z_mean,count,nn_mean,nn_median,nn_p05,nn_p95,pred_spacing,ratio_median_pred`
(`pred_spacing = z_mean * stride / fx`).

Gerçek veri deneyleri için `docs/GERCEK_VERI_PLANI.md`.
