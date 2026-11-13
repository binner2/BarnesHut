#!/usr/bin/env python3
"""Figures for the TUM real-data measurements, drawn from the CSVs only.
   python3 experiments/plot_tum.py <results_dir> <out_prefix> [--inline]
Reads  tum_<seq>_frames.csv, tum_<seq>_info.csv and profiles/*.csv in <results_dir>
and writes <out_prefix>_kde_err.svg, <out_prefix>_meanshift.svg, <out_prefix>_depth_profile.svg.
The SVGs use the report theme classes (s1/s2/s3/s0, l1/l2/l3, lref, grid, tick, lab, lab-b);
a <style> block is embedded so they also render standalone. --inline omits that block
(for embedding into the report, where the theme supplies the classes)."""
import csv, glob, math, os, sys

STYLE = ('<style>.grid{stroke:#d7d9de;stroke-width:1}.tick{fill:#555;font-size:12px}.lab{fill:#222;font-size:13px}'
         '.lab-b{fill:#222;font-size:13px;font-weight:700}.s1{fill:#1f4e9c}.s2{fill:#d9742b}.s3{fill:#2e8b57}.s0{fill:#7a7f88}'
         '.l1{stroke:#1f4e9c;stroke-width:2;fill:none}.l2{stroke:#d9742b;stroke-width:2;fill:none}.l3{stroke:#2e8b57;stroke-width:2;fill:none}'
         '.lref{stroke:#7a7f88;stroke-width:2;fill:none;stroke-dasharray:6 4}.dash{stroke-dasharray:5 4}</style>')
NUM = ('n', 'theta', 'kappa', 'cutoff_sigmas', 'leaf', 'threads', 'build_s', 'time_s', 'far_evals', 'near_evals',
       'max_rel_err', 'mean_rel_err', 'rms_rel_err')

def read_frames(path):
    rows = list(csv.DictReader(open(path)))
    for r in rows:
        for k in NUM:
            if r.get(k, '') != '': r[k] = float(r[k])
        ex = {}
        for kv in (r.get('extra') or '').split(';'):
            if '=' in kv:
                a, b = kv.split('=', 1)
                try: ex[a] = float(b)
                except ValueError: ex[a] = b
        r['ex'] = ex
        r['frame'] = str(ex.get('frame', ''))
    return rows

def sel(rows, **kw):
    out = []
    for r in rows:
        ok = True
        for k, v in kw.items():
            if k.startswith('ex_'):
                x = r['ex'].get(k[3:])
                if x is None or (isinstance(v, float) and abs(float(x) - v) > 1e-9) or (not isinstance(v, float) and x != v): ok = False
            elif isinstance(v, float):
                if abs(r[k] - v) > 1e-9: ok = False
            elif r[k] != v: ok = False
        if ok: out.append(r)
    return out

def frames_of(rows):
    seen = []
    for r in rows:
        if r['frame'] not in seen: seen.append(r['frame'])
    return seen

def fmt(x, d=2): return f'{x:.{d}f}'.replace('.', ',')

def head(W, H, title):
    return [f'<svg class="chart" role="img" viewBox="0 0 {W} {H}" style="min-width:{int(W*0.92)+4}px" xmlns="http://www.w3.org/2000/svg"><title>{title}</title>']

def log_axis_y(o, Y, y0, y1, ml, W, mr, as_pct=True):
    for e in range(int(math.floor(y0)), int(math.ceil(y1)) + 1):
        for m in (1, 2, 5):
            v = m * 10 ** e
            if y0 <= math.log10(v) <= y1:
                o.append(f'<line class="grid" x1="{ml}" y1="{Y(v):.1f}" x2="{W-mr}" y2="{Y(v):.1f}"/>')
                if m == 1:
                    lbl = (fmt(100*v, 2 if v < 1e-3 else (1 if v < 1e-2 else 0)) + ' %') if as_pct else fmt(v, 3)
                    o.append(f'<text class="tick" x="{ml-6}" y="{Y(v)+4:.1f}" text-anchor="end">{lbl}</text>')

def kde_err(rows, out, inline):
    """Per-frame RMS relative error of the four BH settings, one panel per bandwidth."""
    frames = frames_of(rows)
    hs = sorted({r['ex']['h'] for r in sel(rows, algorithm='kde')})
    W, H = 760, 400; ml, mr, mt, mb = 70, 16, 28, 48
    pw = (W - ml - mr - 24) / len(hs)
    configs = [(0.5, 0.0, 'l1', 's1', ''), (0.7, 0.0, 'l1', 's1', 'dash'), (0.5, 1.0, 'l2', 's2', ''), (0.7, 1.0, 'l2', 's2', 'dash')]
    allv = [max(r['rms_rel_err'], 1e-6) for r in sel(rows, algorithm='kde')]
    y0, y1 = math.log10(min(allv)) - 0.2, math.log10(max(allv)) + 0.2
    o = head(W, H, 'TUM freiburg1_xyz: kare başına Barnes-Hut KDE göreli RMS hatası, θ ve κ ayarları, iki bant genişliği')
    if not inline: o.insert(1, STYLE)
    def Y(v): return H - mb - (math.log10(v) - y0) / (y1 - y0) * (H - mt - mb)
    for pi, h in enumerate(hs):
        x_l = ml + pi * (pw + 24)
        def X(i, x_l=x_l): return x_l + (i + 0.5) / len(frames) * pw
        log_axis_y(o, Y, y0, y1, x_l, x_l + pw + mr, mr) if pi == 0 else None
        for e in range(int(math.floor(y0)), int(math.ceil(y1)) + 1):
            v = 10 ** e
            if y0 <= e <= y1: o.append(f'<line class="grid" x1="{x_l}" y1="{Y(v):.1f}" x2="{x_l+pw:.1f}" y2="{Y(v):.1f}"/>')
        for i in range(0, len(frames), 5):
            o.append(f'<text class="tick" x="{X(i):.1f}" y="{H-mb+16}" text-anchor="middle">{i+1}</text>')
        o.append(f'<text class="lab-b" x="{x_l+pw/2:.1f}" y="{mt-10}" text-anchor="middle">h = {fmt(100*h,0)} cm</text>')
        for th, ka, lc, sc, dash in configs:
            pts = []
            for i, f in enumerate(frames):
                r = sel(rows, algorithm='kde', theta=th, kappa=ka, ex_h=h, frame=f)
                if r: pts.append((i, max(r[0]['rms_rel_err'], 1e-6), r[0]['time_s']))
            if not pts: continue
            o.append(f'<polyline class="{lc} {dash}" points="' + ' '.join(f'{X(i):.1f},{Y(v):.1f}' for i, v, _ in pts) + '"/>')
            for i, v, t in pts:
                o.append(f'<circle class="{sc}" cx="{X(i):.1f}" cy="{Y(v):.1f}" r="3.5"><title>kare {i+1}, θ={fmt(th,1)}, κ={fmt(ka,0)}, h={fmt(100*h,0)} cm: RMS {fmt(100*v,3)} %, {fmt(t,3)} s</title></circle>')
    o.append(f'<text class="lab-b" x="{(ml+W-mr)/2:.0f}" y="{H-8}" text-anchor="middle">kare (zaman sırasında, 1–{len(frames)})</text>')
    o.append(f'<text class="lab-b" transform="translate(14,{(mt+H-mb)/2:.0f}) rotate(-90)" text-anchor="middle">göreli RMS hata, log</text>')
    o.append('</svg>')
    open(out, 'w').write('\n'.join(o))

def meanshift(rows, out, inline):
    """Per-frame label agreement with the exact run for theta = 0.5, 0.7, 1.0, and the cluster counts."""
    frames = frames_of(rows)
    W, H = 760, 380; ml, mr, mt, mb = 64, 16, 20, 48
    ths = [(0.5, 'l1', 's1'), (0.7, 'l3', 's3'), (1.0, 'l2', 's2')]
    ag = [r['ex']['agreement'] for r in sel(rows, algorithm='meanshift') if r['theta'] > 0]
    y0 = max(0.0, math.floor(min(ag) * 20) / 20 - 0.05); y1 = 1.0
    def X(i): return ml + (i + 0.5) / len(frames) * (W - ml - mr)
    def Y(v): return H - mb - (v - y0) / (y1 - y0) * (H - mt - mb)
    o = head(W, H, 'TUM freiburg1_xyz: kare başına mean shift etiket uyumu, θ = 0,5, 0,7 ve 1,0')
    if not inline: o.insert(1, STYLE)
    step = 0.05 if (y1 - y0) <= 0.5 else 0.1
    v = y0
    while v <= y1 + 1e-9:
        o.append(f'<line class="grid" x1="{ml}" y1="{Y(v):.1f}" x2="{W-mr}" y2="{Y(v):.1f}"/><text class="tick" x="{ml-6}" y="{Y(v)+4:.1f}" text-anchor="end">{fmt(100*v,0)} %</text>')
        v += step
    for i in range(0, len(frames), 5):
        o.append(f'<text class="tick" x="{X(i):.1f}" y="{H-mb+16}" text-anchor="middle">{i+1}</text>')
    for th, lc, sc in ths:
        pts = []
        for i, f in enumerate(frames):
            r = sel(rows, algorithm='meanshift', theta=th, frame=f)
            if r: pts.append((i, r[0]['ex']['agreement'], int(r[0]['ex']['clusters']), r[0]['time_s']))
        if not pts: continue
        o.append(f'<polyline class="{lc}" points="' + ' '.join(f'{X(i):.1f},{Y(a):.1f}' for i, a, _, _ in pts) + '"/>')
        for i, a, c, t in pts:
            ex = sel(rows, algorithm='meanshift', theta=0.0, frame=frames[i])
            c0 = int(ex[0]['ex']['clusters']) if ex else 0
            o.append(f'<circle class="{sc}" cx="{X(i):.1f}" cy="{Y(a):.1f}" r="4"><title>kare {i+1}, θ={fmt(th,1)}: uyum {fmt(100*a,1)} %, {c} küme (kesin {c0}), {fmt(t,2)} s</title></circle>')
    o.append(f'<text class="lab-b" x="{(ml+W-mr)/2:.0f}" y="{H-8}" text-anchor="middle">kare (zaman sırasında, 1–{len(frames)})</text>')
    o.append(f'<text class="lab-b" transform="translate(14,{(mt+H-mb)/2:.0f}) rotate(-90)" text-anchor="middle">etiket uyumu (kesin çalıştırmaya göre)</text>')
    o.append('</svg>')
    open(out, 'w').write('\n'.join(o))

def read_profiles(results_dir):
    rows = []
    for p in sorted(glob.glob(os.path.join(results_dir, 'profiles', '*.csv'))):
        for r in csv.DictReader(open(p)):
            for k in ('z_lo', 'z_hi', 'z_mean', 'count', 'nn_mean', 'nn_median', 'nn_p05', 'nn_p95', 'pred_spacing', 'ratio_median_pred'):
                r[k] = float(r[k])
            rows.append(r)
    return rows

def fit_slope(prof, min_count=200):
    """Weighted least squares of nn_median = a * z through the origin (weights = bin counts)."""
    num = sum(r['count'] * r['z_mean'] * r['nn_median'] for r in prof if r['count'] >= min_count)
    den = sum(r['count'] * r['z_mean'] ** 2 for r in prof if r['count'] >= min_count)
    return num / den if den else float('nan')

def depth_profile(prof, out, inline, f=525.0, stride=1.0):
    """Median NN spacing per 25 cm depth bin, all frames, against the pinhole prediction z*stride/f."""
    pts = [r for r in prof if r['count'] >= 200]
    W, H = 760, 400; ml, mr, mt, mb = 70, 16, 20, 48
    zmax = max(r['z_hi'] for r in pts); zmin = min(r['z_lo'] for r in pts)
    ymax = max(max(r['nn_p95'] for r in pts), zmax * stride / f) * 1000 * 1.05
    def X(z): return ml + (z - zmin) / (zmax - zmin) * (W - ml - mr)
    def Y(mm): return H - mb - mm / ymax * (H - mt - mb)
    nfr = len(set(r['tag'] for r in prof))
    o = head(W, H, f'TUM freiburg1_xyz: en yakın komşu uzaklığı ile derinlik, {nfr} karenin 25 cm derinlik kutuları, piksel aralığı tahmini z/f')
    if not inline: o.insert(1, STYLE)
    z = math.ceil(zmin * 2) / 2
    while z <= zmax + 1e-9:
        o.append(f'<line class="grid" x1="{X(z):.1f}" y1="{mt}" x2="{X(z):.1f}" y2="{H-mb}"/><text class="tick" x="{X(z):.1f}" y="{H-mb+16}" text-anchor="middle">{fmt(z,1)}</text>')
        z += 0.5
    mm = 0
    while mm <= ymax:
        o.append(f'<line class="grid" x1="{ml}" y1="{Y(mm):.1f}" x2="{W-mr}" y2="{Y(mm):.1f}"/><text class="tick" x="{ml-6}" y="{Y(mm)+4:.1f}" text-anchor="end">{fmt(mm,0)}</text>')
        mm += 1
    o.append(f'<line class="lref" x1="{X(zmin):.1f}" y1="{Y(zmin*stride/f*1000):.1f}" x2="{X(zmax):.1f}" y2="{Y(zmax*stride/f*1000):.1f}"/>')
    o.append(f'<text class="lab" x="{X(zmax)-6:.1f}" y="{Y(zmax*stride/f*1000)-8:.1f}" text-anchor="end">z / f, f = {fmt(f,0)} px</text>')
    for r in pts:
        o.append(f'<line class="l1" style="stroke-width:1;opacity:.35" x1="{X(r["z_mean"]):.1f}" y1="{Y(r["nn_p05"]*1000):.1f}" x2="{X(r["z_mean"]):.1f}" y2="{Y(r["nn_p95"]*1000):.1f}"/>')
        o.append(f'<circle class="s1" cx="{X(r["z_mean"]):.1f}" cy="{Y(r["nn_median"]*1000):.1f}" r="3.5" style="opacity:.8"><title>kare {r["tag"]}, z {fmt(r["z_lo"],2)}–{fmt(r["z_hi"],2)} m, {int(r["count"])} nokta: NN medyan {fmt(r["nn_median"]*1000,2)} mm, tahmin {fmt(r["pred_spacing"]*1000,2)} mm, oran {fmt(r["ratio_median_pred"],3)}</title></circle>')
    a = fit_slope(prof)
    o.append(f'<text class="lab" x="{ml+8}" y="{mt+14}">ağırlıklı uydurma: NN = z / {fmt(1/a,1)} px (çizgi: z / {fmt(f,0)})</text>')
    o.append(f'<text class="lab-b" x="{(ml+W-mr)/2:.0f}" y="{H-8}" text-anchor="middle">derinlik z [m]</text>')
    o.append(f'<text class="lab-b" transform="translate(14,{(mt+H-mb)/2:.0f}) rotate(-90)" text-anchor="middle">NN uzaklığı [mm]: medyan ve 5–95 persentil</text>')
    o.append('</svg>')
    open(out, 'w').write('\n'.join(o))

if __name__ == '__main__':
    d, pre = sys.argv[1], sys.argv[2]
    inline = '--inline' in sys.argv
    fr = glob.glob(os.path.join(d, 'tum_*_frames.csv'))[0]
    rows = read_frames(fr)
    kde_err(rows, pre + '_kde_err.svg', inline)
    meanshift(rows, pre + '_meanshift.svg', inline)
    prof = read_profiles(d)
    if prof: depth_profile(prof, pre + '_depth_profile.svg', inline)
    print('frames:', len(frames_of(rows)), 'rows:', len(rows), 'profile bins:', len(prof), 'fit 1/slope:', fit_slope(prof) and 1 / fit_slope(prof))
