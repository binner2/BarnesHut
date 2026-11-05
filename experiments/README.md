# Deneyler

- `run_benchmarks.sh <pointcloud_tool> <çıktı_dizini>`: raporun bütün sentetik
  ölçümlerini (N × θ × κ taraması, yaprak kapasitesi, uçtan uca demo) yeniden üretir
  ve makine bilgisini `machine.txt`'ye yazar.
- `plot_benchmarks.py <sonuç_dizini> <önek>`: CSV'den Pareto ve ölçeklenme SVG'leri.
- `render_views.py <girdi.ply> <çıktı.webp> [splat]`: renkli PLY'yi sentetik kamera
  iç parametreleriyle görüntüye yeniden izdüşürür (rapordaki 4 panel).
- `results/2026-11-04/`: raporda kullanılan ham ölçümler (4 iş parçacığı, GCC 13.3,
  -O3 -march=native, bulut konteyneri; mutlak süreler donanıma bağlıdır).

CSV şeması: `algorithm,n,theta,kappa,cutoff_sigmas,leaf,threads,build_s,time_s,
far_evals,near_evals,max_rel_err,mean_rel_err,rms_rel_err,extra`.
`extra` alanı `anahtar=değer;...` biçimindedir (mean shift: `clusters`, `agreement`;
WLOP: `nn_cov`, `m`).

Gerçek veri deneyleri için `docs/GERCEK_VERI_PLANI.md`.
