Rapor kaynağı. Yeniden üretim:

    python3 docs/report_src/assemble_report.py docs/report_src/report_template.src.html \
        docs/derinlik_kamerasi_barnes_hut.src.html experiments/results/2026-11-04 <img_dizini>
    node <premium-html>/scripts/build.mjs docs/derinlik_kamerasi_barnes_hut.src.html \
        docs/derinlik_kamerasi_barnes_hut_coworks.html

`<img_dizini>` içindeki webp görüntüler `experiments/render_views.py` ile demo
PLY dosyalarından üretilir (scene_labels, density_full, meanshift_ds, wlop_ds).
