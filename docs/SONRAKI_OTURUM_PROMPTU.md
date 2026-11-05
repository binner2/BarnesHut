# Sonraki Oturum İçin Prompt

Yeni bir Claude Code oturumu açarken aşağıdaki metni olduğu gibi yapıştır.
İş büyük olduğu için dört ayrı oturum önerilir; her oturum için ayrı prompt
verilmiştir. Hepsi bu depodaki `docs/YOL_HARITASI.md` ve
`docs/GERCEK_VERI_PLANI.md` dosyalarını tek doğruluk kaynağı sayar.

---

## Ortak giriş (her prompt'un başına ekle)

```
Depo: binner2/BarnesHut. Önce `git fetch origin claude/sweet-faraday-fyldry` yap ve
bu dalın üzerinden çalış (henüz main'e birleşmedi). İlk iş olarak şu dosyaları oku:
docs/YOL_HARITASI.md, docs/GERCEK_VERI_PLANI.md, pointcloud/README.md,
experiments/README.md. Derle ve testleri çalıştır:
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j && ctest --test-dir build

Kurallar:
- Türkçe yaz, teknik terimler İngilizce, akademik ton.
- Hiçbir sayıyı tahmin etme; her sayı bir CSV'den ya da logdan gelsin. Ölçemediğini
  "ölçülmedi" diye yaz.
- Doğrulayamadığın künye/DOI/URL'yi [DOĞRULANMADI] olarak işaretle; DOI uydurma.
- Bütün sonuçlar experiments/results/<tarih>/ altına, experiments/README.md'deki CSV
  şemasıyla; her CSV'nin yanında machine.txt (CPU, iş parçacığı, derleyici, commit).
- Şekiller yalnızca betikle (experiments/plot_benchmarks.py veya türevi) üretilir.
- Her anlamlı adımda commit at; oturum sonunda aynı dala push et. PR açma.
- Bittiğinde docs/YOL_HARITASI.md'deki durum tablosunu güncelle.
```

---

## Oturum 1 — Gerçek veri: TUM RGB-D (1. öncelik)

```
Görev: docs/GERCEK_VERI_PLANI.md bölüm 0, 1.1 ve 3'ü uygula.

1. TUM RGB-D freiburg1_xyz dizisini indir
   (https://cvg.cit.tum.de/rgbd/dataset/freiburg1/rgbd_dataset_freiburg1_xyz.tgz).
   Disk kısıtlıysa yalnızca depth/ klasörünü tut. İndirme engellenirse dur ve bana
   hangi dosyayı nereye koymam gerektiğini söyle.
2. tools/tum_depth_to_ply.sh ile derinlik PNG'lerini PLY'ye çevir (ölçek 5000,
   fx=fy=525, cx=319.5, cy=239.5). İlk kareyi `pointcloud_tool info` ile incele ve
   nokta sayısı, sınır kutusu, NN uzaklık istatistiklerini kaydet.
3. 20 kare seç (eşit aralıklı). Her kare için şu ölçümleri al ve
   experiments/results/<tarih>/tum_freiburg1_xyz_frames.csv'ye yaz:
   octree inşa süresi; KDE kesin (theta=0, kesme yok) vs theta∈{0.5,0.7} ×
   kappa∈{0,1} için süre ve göreli RMS/max hata (h = 2 cm ve 5 cm); mean shift
   (h=8 cm) theta∈{0,0.5,0.7} süre ve etiket uyumu; WLOP (M = N/10) theta∈{0,0.5,0.7}
   süre ve NN CoV. 5 tekrar, medyan raporla.
4. Sentetik sahnedeki bulgu gerçek veride tutuyor mu? Özellikle: kappa=1 geometrik
   ölçüte göre hatayı kaç kat düşürüyor; theta=1.0 mean shift'i gerçek veride de
   bozuyor mu; yoğunluğun 1/z^2 ile değişimi (eq. rho_A) ölçülen NN uzaklığıyla
   uyumlu mu? Her soruya CSV'ye dayanan tek paragraf yanıt ver.
5. Sonuçları docs/report_src/report_template.src.html'e yeni bir "Gerçek veri"
   bölümü olarak ekle, raporu yeniden build et (docs/report_src/README.md) ve
   "Sınırlılıklar" bölümünü güncelle.

Teslim: CSV'ler, güncel rapor, YOL_HARITASI durum tablosu, push.
```

---

## Oturum 2 — Taban çizgileri (2. öncelik)

```
Görev: docs/GERCEK_VERI_PLANI.md bölüm 2'deki taban çizgilerini kur ve Oturum 1'in
20 karesinde aynı parametrelerle çalıştır.

- pip ile open3d kur (sürümü kaydet). Open3D KD-tree yarıçap aramasıyla aynı h için
  Gauss KDE'yi Python'da hesapla ve süreyi ölç; voxel_down_sample,
  farthest_point_down_sample ile seyreltme; sample_points_poisson_disk için ağ
  gerekiyorsa bunu belirt ve atla.
- CGAL kurulabiliyorsa wlop_simplify_and_regularize_point_set ile WLOP; aynı h, mu,
  yineleme. Kurulamıyorsa nedenini yaz ve atla.
- IFGT/FIGTree kaynak koddan derlenebiliyorsa aynı h ve hata hedefiyle KDE;
  derlenemiyorsa yaz ve atla.
- Bütün sonuçlar experiments/results/<tarih>/baselines.csv'ye (yöntem, kütüphane,
  sürüm, kare, h, süre, hata, CoV).
- experiments/plot_benchmarks.py'ye "bizim yöntem vs taban çizgileri" Pareto
  grafiği ekle ve rapora işle.

Teslim: baselines.csv, güncel rapor, push. Hiçbir taban çizgisini "tahminen" yazma.
```

---

## Oturum 3 — Eksen A teorisi ve kod eksikleri (3–4. öncelik)

```
Görev A (teori): docs/derinlik_kamerasi_barnes_hut_coworks.html bölüm 4.3'teki
ikinci mertebe hata kestirimini (eq-err) üçüncü mertebe kalanı da içeren kesin bir
üst sınıra dönüştür. Türetmeyi adım atlamadan yaz, Gauss çekirdeği için Taylor
kalanını Lagrange biçiminde sınırla. Sonra bu sınırı sentetik benchmark
(experiments/results/2026-11-04/bench_main.csv) ve varsa TUM sonuçlarıyla
karşılaştır: ölçülen en büyük hata sınırın altında mı, sınır ne kadar gevşek?
theta, kappa ve yaprak kapasitesi için h/nokta aralığı oranına bağlı bir seçim
kuralı öner ve bu kuralı pointcloud_tool'a `--auto` seçeneği olarak ekle.

Görev B (kod): docs/GERCEK_VERI_PLANI.md bölüm 4'teki maddeleri sırayla uygula:
değişken bant genişliği h(z); mean shift için adım normuna göre erken durdurma;
konum–normal 6-D mean shift (düğümlerde normal ortalaması, açılma ölçütüne normal
dağılımı); WLOP çekim kesmesi 1.5h; PCD binary_compressed (LZF) okuma. Her madde
için tests/test_pointcloud.cpp'ye test ekle; ctest geçmeden commit atma.

Teslim: türetme raporun yeni bir bölümü olarak (KaTeX), kod + testler, push.
```

---

## Oturum 4 — Makale taslağı (5. öncelik; 1–3 bittikten sonra)

```
Görev: IEEE konferans biçiminde (iki sütun, LaTeX, IEEEtran) makale taslağı yaz.
Konu: Eksen A (çekirdek-bilinçli Barnes-Hut açılma ölçütü) + Eksen B (BH-WLOP ile
derinlik kamerası örnekleme yanlılığının giderilmesi). Kaynaklar:
docs/derinlik_kamerasi_barnes_hut_coworks.html (yöntem, türetme),
experiments/results/ altındaki bütün CSV'ler (deneyler), docs/report_src (şekiller).

Adımlar:
1. Önce ieee-cs-writing ve literatur-kill-search becerilerini kullan: "kernel-aware
   opening criterion octree Gaussian kernel density", "Barnes-Hut kernel density
   estimation point cloud", "WLOP depth camera sampling bias" sorgularıyla yakın
   rakipleri bul; bulduklarını related work'e koy, bulunamayanı [DOĞRULANMADI]
   işaretle. Katkı iddialarını bu taramaya göre daralt.
2. Yapı: Abstract, Introduction (katkı listesi 3 madde), Related Work, Method
   (octree, BH toplamları, hata sınırı, kappa ölçütü, WLOP), Experiments (sentetik
   + TUM + taban çizgileri; her şekil CSV'den), Limitations, Conclusion.
3. Bütün sayılar CSV'den; şekiller experiments/plot_benchmarks.py'den PDF olarak.
4. paper-contribution-gate becerisiyle katkı iddialarını denetle; overclaiming varsa
   düzelt.
5. paper/ dizinine LaTeX kaynakları, derlenmiş PDF ve derleme komutu.

Teslim: paper/main.pdf, kill-search günlüğü, push.
```

---

## Pratik notlar

- Oturumlar bağımsız dallar değil; hepsi `claude/sweet-faraday-fyldry` üzerinde
  sırayla çalışsın ki sonuçlar birikerek rapora girsin. İstersen önce bu dalı main'e
  birleştir ve prompt'lardaki dal adını kaldır.
- TUM dizisi ~0,5 GB; bulut konteynerinde disk kotası dolarsa yalnızca depth/
  klasörünü tut, stride 2 kullan.
- Open3D'nin CUDA'sız pip paketi yeterlidir; CGAL ve FIGTree kurulamazsa atlanır ve
  raporda "ölçülmedi" denir.
- Kendi RealSense/Azure Kinect kaydın varsa Oturum 1'e ekle:
  "data/own/<sahne>/ altında 16-bit PNG + intrinsics.txt var; aynı ölçümleri bu
  karelerde de al." (tools/realsense_capture.py biçimi.)
