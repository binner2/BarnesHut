# Barnes-Hut N-Body Simülasyonu

[![C++](https://img.shields.io/badge/Language-C++-blue.svg)](https://isocpp.org/)
[![Algorithm](https://img.shields.io/badge/Algorithm-Barnes--Hut-green.svg)](https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation)

## 📋 İçindekiler

- [Proje Hakkında](#proje-hakkında)
- [Barnes-Hut Algoritması](#barnes-hut-algoritması)
- [Versiyon Geçmişi](#versiyon-geçmişi)
- [Proje Yapısı](#proje-yapısı)
- [Kurulum ve Derleme](#kurulum-ve-derleme)
- [Kullanım](#kullanım)
- [Görselleştirme İmkanları](#görselleştirme-imkanları)
- [Parametreler](#parametreler)
- [Teknik Detaylar](#teknik-detaylar)

---

## 🎯 Proje Hakkında

Bu proje, **N-Body problemi**ni çözmek için **Barnes-Hut algoritması** kullanan bir gravitasyonel/elektromanyetik simülasyon uygulamasıdır. N-Body problemi, birbirleriyle etkileşen N adet cismin (parçacık, yıldız, galaksi vb.) zaman içindeki hareketlerini hesaplamayı amaçlar.

### Temel Özellikler

- **Barnes-Hut Octree** yapısı ile hızlandırılmış kuvvet hesaplaması
- **Leapfrog entegrasyon** yöntemi ile zaman adımlaması
- **3 boyutlu** uzayda parçacık simülasyonu
- **Coulomb** veya **gravitasyonel** kuvvet hesaplaması
- Performans ölçümleme ve zaman takibi
- Gnuplot ile uyumlu çıktı formatı
- Özelleştirilebilir TETA ve LEAF parametreleri

---

## 🌳 Barnes-Hut Algoritması

Barnes-Hut algoritması, N-Body simülasyonlarında hesaplama karmaşıklığını **O(N²)** den **O(N log N)** seviyesine indiren bir hiyerarşik ağaç tabanlı yaklaşımdır.

### Algoritma Prensibi

1. **Ağaç Oluşturma**: 3D uzay, recursive olarak 8 alt kübe (octree) bölünür
2. **Kütle Merkezi Hesaplama** (Upward Pass): Her düğüm için toplam kütle ve kütle merkezi hesaplanır
3. **Kuvvet Hesaplama**: Uzak parçacık grupları tek bir kütle merkezi olarak yaklaşıklanır
4. **TETA Kriteri**: Bir hücrenin ne kadar uzakta olduğunu belirler
   - `s/d < TETA` ise hücre tek nokta gibi işlem görür
   - `s/d >= TETA` ise hücre alt parçalarına ayrılır
   - (s: hücre boyutu, d: mesafe)

### Avantajları

- ✅ Hesaplama süresinde dramatik azalma (özellikle büyük N değerlerinde)
- ✅ Bellek kullanımında verimlilik
- ✅ Ayarlanabilir doğruluk-performans dengesi (TETA parametresi)
- ✅ Büyük ölçekli simülasyonları mümkün kılar

---

## 📜 Versiyon Geçmişi

### 🔵 İlk Sürüm (Mart 2017 - Commit: 1150034)

**Yazar**: Burak İNNER
**Tarih**: 8 Mart 2017

İlk implementasyon şu özellikleri içeriyordu:

#### Ana Bileşenler
- **BHtreetest.cpp**: Ana program ve simülasyon döngüsü
- **tree.cpp/h**: Barnes-Hut ağaç yapısı ve algoritma implementasyonu
  - Octree oluşturma
  - Upward pass (kütle merkezi hesaplama)
  - Downward pass (kuvvet hesaplama)
  - Parçacık yükleme ve ağaç yönetimi
- **particle.cpp/h**: Parçacık sınıfı
  - Pozisyon, hız, kuvvet, kütle bilgileri
  - Leapfrog entegrasyonu
- **vektor.h**: 3D vektör matematik kütüphanesi
  - Vektör operatörleri (+, -, *, /)
  - Dot product, mesafe hesaplama
- **file.cpp/h**: Dosya I/O işlemleri
  - Parçacık verilerini okuma
  - Sonuçları yazma (pozisyon ve kuvvet)

#### Temel Özellikler
- Fast Multipole Method (FMM) yaklaşımı
- Coulomb kuvveti hesaplaması
- Zaman ölçümleme (#define showtime)
- Leapfrog entegrasyon şeması
- Parametrik TETA ve particleLEAF değerleri

#### Geliştirme Notları (Kod İçi Yorumlardan)
- **07-02-2005**: Kod gözden geçirildi
- **01-03-2005**:
  - TETA ve particleLEAF global değişkenler yapıldı
  - LEAF yapısı STL::vector kullanmaya başladı
  - Nodelist STL::vector'e dönüştürüldü
- **26-08-2004**: TETA parametresi static'ten fonksiyon argümanına çevrildi

### 🟢 Güncel Sürüm (Nisan 2018 - Commit: ef3a95b)

**Yazar**: yapbenzet
**Tarih**: 11 Nisan 2018

#### Yeni Özellikler
- **generate_data.cpp**: Rastgele test verisi üreteci eklendi
  - Kullanıcı tanımlı parçacık sayısı
  - Rastgele kütle, pozisyon ve hız değerleri
  - Simülasyon parametrelerini tanımlama
  - Standart veri formatında dosya oluşturma

#### Farklar ve İyileştirmeler

| Özellik | İlk Sürüm | Güncel Sürüm |
|---------|-----------|--------------|
| Veri Üretimi | Manuel dosya oluşturma | Otomatik rastgele veri üreteci |
| Test Kolaylığı | Harici araçlar gerekli | Tek komutla test verisi |
| Esneklik | Sabit veri setleri | Parametrik veri üretimi |
| Kullanım | Veri dosyası hazırlama gerekli | generate_data ile hızlı başlangıç |

---

## 📁 Proje Yapısı

```
BarnesHut/
│
├── BHtreetest.cpp          # Ana program (main fonksiyonu)
├── tree.cpp / tree.h       # Barnes-Hut ağaç implementasyonu
├── particle.cpp / particle.h   # Parçacık sınıfı tanımı
├── vektor.h                # 3D vektör matematik kütüphanesi
├── file.cpp / file.h       # Dosya I/O işlemleri
├── stdinc.cpp / stdinc.h   # Standart include ve sabitler
├── generate_data.cpp       # Rastgele veri üreteci
└── ncb                     # Derlenmiş binary (executable)
```

### Dosya Açıklamaları

#### **BHtreetest.cpp**
- Program giriş noktası
- Komut satırı argümanlarını işler
- Ana simülasyon döngüsü
- Zaman adımlarını yönetir
- Performans metrikleri toplar

#### **tree.cpp/h**
- Octree veri yapısı
- Parçacık yükleme (loadparticles)
- Upward pass (upwardpass): Kütle merkezi hesaplama
- Downward pass (downwardpass): Multipole-to-local dönüşümü
- Kuvvet hesaplama (calculateforces)
- Ağaç oluşturma ve temizleme

#### **particle.cpp/h**
- Particle sınıfı: Kütle, pozisyon, hız, kuvvet
- Node struct: Octree düğüm yapısı
- Leapfrog entegrasyon (calcVelPos)
- Parçacık bilgilerini görüntüleme

#### **vektor.h**
- 3D vektör sınıfı
- Matematiksel operatörler
- Dot product, mesafe fonksiyonları
- Inline optimizasyonlar

#### **file.cpp/h**
- readfilePart: Parçacık verilerini dosyadan okuma
- putsnapshotfilePos: Pozisyon snapshot'ları
- putsnapshotfileForce: Kuvvet snapshot'ları
- Gnuplot uyumlu format

#### **generate_data.cpp**
- İnteraktif veri üretici
- Rastgele kütle ataması (5000-15000 arası)
- Rastgele pozisyon (0-10 arası)
- Rastgele hız (0-100 arası)

---

## 🔧 Kurulum ve Derleme

### Gereksinimler

- **C++ Derleyici**: g++ (GCC 4.x veya üzeri)
- **Standart Kütüphaneler**: STL (vector, iostream, fstream)
- **İşletim Sistemi**: Linux, macOS veya Windows (MinGW)

### Derleme Adımları

#### 1. Ana Programı Derleme

```bash
# Tüm kaynak dosyaları birlikte derle
g++ -o bh_simulation BHtreetest.cpp tree.cpp particle.cpp file.cpp stdinc.cpp -lm -O3

# veya optimize edilmiş versiyon
g++ -o bh_simulation BHtreetest.cpp tree.cpp particle.cpp file.cpp stdinc.cpp \
    -lm -O3 -march=native -funroll-loops
```

#### 2. Veri Üreteciyi Derleme

```bash
g++ -o generate_data generate_data.cpp -O2
```

### Derleme Parametreleri

- `-O3`: Maksimum optimizasyon
- `-march=native`: İşlemciye özel optimizasyon
- `-funroll-loops`: Döngü açma optimizasyonu
- `-lm`: Matematik kütüphanesi

---

## 🚀 Kullanım

### Adım 1: Test Verisi Oluşturma

İki yöntem mevcuttur:

#### Yöntem A: Otomatik Veri Üreteci (Önerilen)

```bash
./generate_data
```

Program sizden şu bilgileri isteyecektir:
- Dosya adı (örn: `test_data.dat`)
- Parçacık sayısı (örn: `1000`)
- Başlangıç zamanı (örn: `0.0`)
- Bitiş zamanı (örn: `10.0`)
- Zaman adımı (örn: `0.01`)

Veri otomatik olarak oluşturulur.

#### Yöntem B: Manuel Veri Dosyası

`input.dat` dosyası şu formatta olmalıdır:

```
N                           # Parçacık sayısı
t_start                     # Başlangıç zamanı
t_end                       # Bitiş zamanı
dt                          # Zaman adımı
mass[1] x[1] y[1] z[1] vx[1] vy[1] vz[1]
mass[2] x[2] y[2] z[2] vx[2] vy[2] vz[2]
...
mass[N] x[N] y[N] z[N] vx[N] vy[N] vz[N]
```

### Adım 2: Simülasyonu Çalıştırma

```bash
./bh_simulation <veri_dosyası> <TETA> <particleLEAF>
```

#### Örnek Kullanımlar

```bash
# Küçük sistem, yüksek doğruluk
./bh_simulation test_data.dat 0.5 1

# Orta ölçek, dengeli parametreler
./bh_simulation test_data.dat 0.7 8

# Büyük sistem, hızlı hesaplama
./bh_simulation test_data.dat 1.0 16
```

### Adım 3: Çıktıları İnceleme

Program şu bilgileri görüntüler:

```
Starting a BARNES HUT algorithm with Lepfrog integration for a
1000-body system, from time t=0 to t=10 with time step=0.01
and number of max. particle in a LEAF 8 TETA=0.7

load particles takes 0.023 seconds
upward pass (without load particles) takes 0.015 seconds
calculate forces take (without load particles) 0.087 seconds
all takes 0.125 seconds
number of step = 1
```

---

## 📊 Görselleştirme İmkanları

Proje, sonuçları Gnuplot ile görselleştirmek için uygun formatta çıktılar üretir.

### 1. Çıktı Dosyaları

Program iki tür snapshot dosyası oluşturur:

#### Pozisyon Çıktısı
```
putsnapshotfilePos(N, partList, message, TETA);
```
Format: `x y z` (her satırda bir parçacık)

#### Kuvvet Çıktısı
```
putsnapshotfileForce(N, partList, timetakes);
```
Format: `fx fy fz` (her satırda bir parçacığın kuvveti)

### 2. Gnuplot ile Görselleştirme

#### 3D Parçacık Dağılımı

```gnuplot
# Gnuplot betiği - plot_positions.gnu
set terminal png size 1024,768
set output 'particle_distribution.png'

splot 'positions.dat' using 1:2:3 with points pt 7 ps 0.5 title 'Particles'
```

Çalıştırma:
```bash
gnuplot plot_positions.gnu
```

#### 2D Projeksiyonlar

```gnuplot
# XY düzlemi
plot 'positions.dat' using 1:2 with points

# XZ düzlemi
plot 'positions.dat' using 1:3 with points

# YZ düzlemi
plot 'positions.dat' using 2:3 with points
```

#### Animasyon Oluşturma

```bash
# Her zaman adımında snapshot al
# Gnuplot ile frame'leri birleştir

set terminal gif animate delay 10
set output 'simulation.gif'

do for [i=0:100] {
    splot sprintf('snapshot_%03d.dat', i) using 1:2:3 with points
}
```

#### Kuvvet Vektörleri

```gnuplot
# Kuvvet vektörlerini oklar olarak göster
plot 'data.dat' using 1:2:3:4 with vectors
```

### 3. Veri Analizi İmkanları

#### Enerji Takibi
Kinetik ve potansiyel enerjileri izleyerek sistemin korunumunu kontrol edin.

#### Yörünge Analizi
Belirli parçacıkların yörüngelerini zaman içinde takip edin.

#### Dağılım İstatistikleri
Parçacık dağılımının zaman içindeki evrimini inceleyin.

---

## ⚙️ Parametreler

### TETA (θ) Parametresi

Barnes-Hut algoritmasının doğruluk-performans dengesini belirler.

| TETA Değeri | Doğruluk | Hız | Kullanım Durumu |
|-------------|----------|-----|-----------------|
| 0.0 - 0.5   | Çok yüksek | Yavaş | Kesin hesaplama gerekli |
| 0.5 - 0.7   | Yüksek | Orta | Genel amaçlı (önerilen) |
| 0.7 - 1.0   | Orta | Hızlı | Büyük sistemler |
| > 1.0       | Düşük | Çok hızlı | Kaba yaklaşım |

**Açıklama**:
- `TETA = s/d` (s: hücre boyutu, d: mesafe)
- Küçük TETA → Daha fazla direkt hesaplama → Daha doğru
- Büyük TETA → Daha fazla yaklaşıklama → Daha hızlı

### particleLEAF Parametresi

Bir LEAF düğümünde bulunabilecek maksimum parçacık sayısı.

| particleLEAF | Ağaç Derinliği | Bellek | Kullanım Durumu |
|--------------|----------------|--------|-----------------|
| 1            | Derin | Fazla | Maksimum doğruluk |
| 4-8          | Dengeli | Orta | Genel amaçlı (önerilen) |
| 16-32        | Sığ | Az | Bellek kısıtlı sistemler |

**Açıklama**:
- Küçük değer → Daha derin ağaç → Daha iyi kuvvet yaklaşıklaması
- Büyük değer → Daha sığ ağaç → Daha az bellek kullanımı

### Zaman Adımı (dt)

Simülasyon zaman çözünürlüğü.

```
dt = 0.001    # Çok küçük adım - kararlı ama yavaş
dt = 0.01     # Orta - çoğu durum için uygun
dt = 0.1      # Büyük adım - hızlı ama kararsızlık riski
```

**Kural**: `dt < 0.1 * (karakteristik zaman ölçeği)`

### Önerilen Parametre Kombinasyonları

```bash
# Hassas bilimsel hesaplama
./bh_simulation data.dat 0.5 1

# Genel amaçlı (1K-10K parçacık)
./bh_simulation data.dat 0.7 8

# Büyük ölçekli (>100K parçacık)
./bh_simulation data.dat 0.9 16
```

---

## 🔬 Teknik Detaylar

### Algoritma Akışı

```
1. Başlangıç
   ├── Veri dosyasından parçacık bilgilerini oku
   ├── TETA ve particleLEAF parametrelerini al
   └── Tree nesnesini oluştur

2. Ana Simülasyon Döngüsü (her zaman adımı için)
   ├── loadparticles()
   │   ├── Ağacı temizle
   │   ├── Root'u bul (tüm parçacıkları içeren küp)
   │   └── Her parçacığı ağaca ekle (addparticle)
   │
   ├── upwardpass()
   │   ├── Yapraktan köke doğru git
   │   └── Her düğüm için kütle merkezi hesapla
   │
   ├── calculateforces()
   │   ├── Her parçacık için
   │   │   ├── Ağacı traverse et
   │   │   ├── TETA kriterine göre karar ver
   │   │   │   ├── Yakın ise → direkt hesaplama
   │   │   │   └── Uzak ise → kütle merkezi yaklaşımı
   │   │   └── Toplam kuvveti hesapla
   │   └── Leapfrog entegrasyonu ile pozisyon/hız güncelle
   │
   ├── emptytree()
   │   └── Ağaç yapısını temizle
   │
   └── Çıktıları dosyaya yaz

3. Sonlandırma
   └── Performans metriklerini raporla
```

### Veri Yapıları

#### Node Yapısı
```cpp
struct Node {
    long index;              // Düğüm indeksi
    _type Type;              // EMPTY, NODE veya LEAF
    vektor geocenter;        // Geometrik merkez
    double rsize;            // Küp kenar uzunluğu
    vektor masscenter;       // Kütle merkezi
    double mass;             // Toplam kütle
    vector<particle*> plist; // LEAF ise parçacık listesi
    int Nparticle;           // Parçacık sayısı
    long level;              // Ağaç seviyesi
    Node* child[8];          // 8 alt düğüm (octree)
    Node* parent;            // Ebeveyn düğüm
};
```

#### Particle Sınıfı
```cpp
class particle {
    double data;       // Kütle veya yük
    vektor pos;        // Pozisyon (x,y,z)
    vektor vel;        // Hız (vx,vy,vz)
    vektor force;      // Kuvvet (fx,fy,fz)
    unsigned int Id;   // Parçacık kimliği
    Node* parent;      // Bağlı olduğu düğüm
};
```

### Leapfrog Entegrasyonu

Parçacık pozisyon ve hızları şu şekilde güncellenir:

```
v(t+dt/2) = v(t) + a(t) * dt/2
x(t+dt)   = x(t) + v(t+dt/2) * dt
v(t+dt)   = v(t+dt/2) + a(t+dt) * dt/2
```

Avantajları:
- Simplektik (enerji korunumu)
- İkinci mertebe doğruluk
- Uzun süreli kararlılık

### Kuvvet Hesaplama

#### Direkt Hesaplama (Yakın Parçacıklar)
```cpp
F = G * m1 * m2 / r² * r̂
```

#### Multipole Yaklaşımı (Uzak Gruplar)
```cpp
F ≈ G * m_particle * M_cell / R² * R̂
```
(M_cell: Hücrenin toplam kütlesi, R: Kütle merkezine mesafe)

### Performans Karakteristikleri

| N (Parçacık) | O(N²) Süresi | O(N log N) Süresi | Hızlanma |
|--------------|--------------|-------------------|----------|
| 100          | 0.01 s       | 0.01 s            | 1x       |
| 1,000        | 1 s          | 0.1 s             | 10x      |
| 10,000       | 100 s        | 1.3 s             | 77x      |
| 100,000      | 10,000 s     | 17 s              | 588x     |

### Bellek Kullanımı

- Her Node: ~120 byte
- Her Particle: ~80 byte
- Toplam ağaç düğümleri: yaklaşık `N/particleLEAF * 1.3` (overhead)

**Örnek**: 10,000 parçacık, particleLEAF=8
- Parçacıklar: 10,000 * 80 = 800 KB
- Düğümler: ~1,600 * 120 = 192 KB
- **Toplam**: ~1 MB

---

## 📚 Referanslar

1. **Barnes, J. & Hut, P.** (1986). "A hierarchical O(N log N) force-calculation algorithm". Nature, 324, 446-449.
2. **Aarseth, S. J.** (2003). "Gravitational N-Body Simulations". Cambridge University Press.
3. **Dehnen, W.** (2002). "A Hierarchical O(N) Force Calculation Algorithm". Journal of Computational Physics, 179, 27-42.

---

## 👨‍💻 Yazar ve Tarihçe

### Orijinal Geliştirici
**Burak İNNER**
Email: binner@kocaeli.edu.tr
Tarih: 6-7 Temmuz 2004 (ilk geliştirme)

### Güncellemeler
- **Mart 2017**: Ana implementasyon
- **Nisan 2018**: Veri üreteci eklendi (yapbenzet@kocaeli.edu.tr)

### Kurum
Kocaeli Üniversitesi

---

## 📄 Lisans

Bu proje eğitim ve araştırma amaçlıdır.

---

## 🤝 Katkıda Bulunma

Projeyi geliştirmek isterseniz:

1. Repository'yi fork edin
2. Feature branch oluşturun (`git checkout -b feature/AmazingFeature`)
3. Değişikliklerinizi commit edin (`git commit -m 'Add some AmazingFeature'`)
4. Branch'inizi push edin (`git push origin feature/AmazingFeature`)
5. Pull Request açın

---

## ❓ Sık Sorulan Sorular

### Hangi TETA değerini seçmeliyim?
- Başlangıç için **0.7** önerilir. Daha hassas sonuçlar için 0.5'e düşürün, hız için 1.0'a çıkarın.

### Program çok yavaş çalışıyor, ne yapmalıyım?
- TETA değerini artırın (0.9-1.0)
- particleLEAF değerini artırın (16-32)
- Derleme optimizasyonlarını kullanın (-O3)

### Sonuçlar kararsız görünüyor?
- Zaman adımını (dt) küçültün
- TETA değerini düşürün
- Parçacık başlangıç koşullarını kontrol edin

### Görselleştirme nasıl yapılır?
- Çıktı dosyalarını Gnuplot ile açın (yukarıdaki örneklere bakın)
- Python matplotlib ile de görselleştirme yapılabilir

---

**Son Güncelleme**: Kasım 2025
