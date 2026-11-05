#!/usr/bin/env python3
"""Draw the report figures from benchmark CSVs as standalone SVG files.
   python3 experiments/plot_benchmarks.py <results_dir> <out_prefix>
Produces <out_prefix>_pareto.svg and <out_prefix>_scaling.svg (log-log)."""
import csv, math, sys, os

C = {'s0': '#7a7f88', 's1': '#1f4e9c', 's2': '#d9742b', 's3': '#2e8b57', 'grid': '#d7d9de', 'txt': '#222'}

def read(path):
    rows = list(csv.DictReader(open(path)))
    for r in rows:
        for k in ('n', 'theta', 'kappa', 'cutoff_sigmas', 'time_s', 'rms_rel_err', 'max_rel_err'):
            r[k] = float(r[k])
    return rows

def sel(rows, **kw):
    return [r for r in rows if all(abs(r[k] - v) < 1e-9 if isinstance(v, float) else r[k] == v for k, v in kw.items())]

def axes(W, H, ml, mr, mt, mb, x0, x1, y0, y1):
    X = lambda v: ml + (math.log10(v) - x0) / (x1 - x0) * (W - ml - mr)
    Y = lambda v: H - mb - (math.log10(v) - y0) / (y1 - y0) * (H - mt - mb)
    return X, Y

def pareto(rows, out):
    nmax = max(int(r['n']) for r in rows)
    ex = sel(rows, algorithm='kde_exact', n=float(nmax))[0]
    pts = []
    for th in (0.3, 0.5, 0.7, 1.0):
        for ka, c in ((0.0, 's1'), (1.0, 's2'), (0.5, 's3')):
            r = sel(rows, algorithm='kde', n=float(nmax), theta=th, kappa=ka, cutoff_sigmas=4.0)
            if r: pts.append((r[0]['time_s'], max(r[0]['rms_rel_err'], 1e-6), f'θ={th}' + ('' if ka == 0 else f', κ={ka}'), c))
    cut = sel(rows, algorithm='kde', n=float(nmax), theta=0.0, cutoff_sigmas=4.0)[0]
    pts.append((cut['time_s'], max(cut['rms_rel_err'], 1e-6), 'kesme', 's0'))
    W, H = 720, 420
    xs = [p[0] for p in pts] + [ex['time_s']]; ys = [p[1] for p in pts]
    X, Y = axes(W, H, 64, 20, 16, 52, math.log10(min(xs)) - .15, math.log10(max(xs)) + .15, math.log10(min(ys)) - .3, math.log10(max(ys)) + .3)
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" font-family="sans-serif" font-size="12" fill="{C["txt"]}">',
         f'<rect width="{W}" height="{H}" fill="white"/>']
    for e in range(-3, 3):
        for m in (1, 2, 5):
            v = m * 10 ** e
            if min(xs) * 0.7 <= v <= max(xs) * 1.4:
                o.append(f'<line x1="{X(v):.1f}" y1="16" x2="{X(v):.1f}" y2="{H-52}" stroke="{C["grid"]}"/><text x="{X(v):.1f}" y="{H-36}" text-anchor="middle">{v:g}</text>')
    for e in range(-6, 1):
        v = 10 ** e
        if min(ys) * 0.5 <= v <= max(ys) * 2:
            o.append(f'<line x1="64" y1="{Y(v):.1f}" x2="{W-20}" y2="{Y(v):.1f}" stroke="{C["grid"]}"/><text x="58" y="{Y(v)+4:.1f}" text-anchor="end">{100*v:g} %</text>')
    o.append(f'<line x1="{X(ex["time_s"]):.1f}" y1="16" x2="{X(ex["time_s"]):.1f}" y2="{H-52}" stroke="{C["s0"]}" stroke-dasharray="6 4"/>')
    o.append(f'<text x="{X(ex["time_s"])-6:.1f}" y="30" text-anchor="end">kesin ({ex["time_s"]:.1f} s)</text>')
    for t, e, lbl, c in pts:
        o.append(f'<circle cx="{X(t):.1f}" cy="{Y(e):.1f}" r="6" fill="{C[c]}"><title>{lbl}: {t:.3f} s, RMS {100*e:.3f} %</title></circle>')
        if c in ('s1', 's0'): o.append(f'<text x="{X(t)+9:.1f}" y="{Y(e)-7:.1f}">{lbl}</text>')
    o.append(f'<text x="{W/2:.0f}" y="{H-10}" text-anchor="middle" font-weight="bold">hesap süresi [s], log</text>')
    o.append(f'<text transform="translate(14,{H/2:.0f}) rotate(-90)" text-anchor="middle" font-weight="bold">göreli RMS hata, log</text>')
    o.append(f'<text x="64" y="12">N={nmax}, h=5 cm; mavi κ kapalı, turuncu κ=1, yeşil κ=0,5, gri yalnız kesme</text></svg>')
    open(out, 'w').write('\n'.join(o))

def scaling(rows, out):
    sizes = sorted({int(r['n']) for r in rows})
    series = [('kesin', 's0', dict(algorithm='kde_exact')), ('θ=0,5', 's1', dict(algorithm='kde', theta=0.5, kappa=0.0, cutoff_sigmas=4.0)),
              ('θ=0,5, κ=1', 's2', dict(algorithm='kde', theta=0.5, kappa=1.0, cutoff_sigmas=4.0)), ('θ=0,7', 's3', dict(algorithm='kde', theta=0.7, kappa=0.0, cutoff_sigmas=4.0))]
    data = {n: [(s, sel(rows, n=float(s), **kw)[0]['time_s']) for s in sizes if sel(rows, n=float(s), **kw)] for n, _, kw in series}
    allt = [t for v in data.values() for _, t in v]
    W, H = 720, 400
    X, Y = axes(W, H, 64, 20, 16, 52, math.log10(min(sizes)) - .1, math.log10(max(sizes)) + .1, math.log10(min(allt)) - .2, math.log10(max(allt)) + .2)
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" font-family="sans-serif" font-size="12" fill="{C["txt"]}">', f'<rect width="{W}" height="{H}" fill="white"/>']
    for s in sizes:
        o.append(f'<line x1="{X(s):.1f}" y1="16" x2="{X(s):.1f}" y2="{H-52}" stroke="{C["grid"]}"/><text x="{X(s):.1f}" y="{H-36}" text-anchor="middle">{s//1000}k</text>')
    for e in range(-3, 4):
        v = 10 ** e
        if min(allt) * .5 <= v <= max(allt) * 2:
            o.append(f'<line x1="64" y1="{Y(v):.1f}" x2="{W-20}" y2="{Y(v):.1f}" stroke="{C["grid"]}"/><text x="58" y="{Y(v)+4:.1f}" text-anchor="end">{v:g} s</text>')
    for name, c, _ in series:
        v = data[name]
        o.append(f'<polyline fill="none" stroke="{C[c]}" stroke-width="2" points="' + ' '.join(f'{X(n):.1f},{Y(t):.1f}' for n, t in v) + '"/>')
        for n, t in v: o.append(f'<circle cx="{X(n):.1f}" cy="{Y(t):.1f}" r="5" fill="{C[c]}"><title>{name}, N={n}: {t:.3f} s</title></circle>')
        o.append(f'<text x="{X(v[-1][0])+8:.1f}" y="{Y(v[-1][1])+4:.1f}">{name}</text>')
    o.append(f'<text x="{W/2:.0f}" y="{H-10}" text-anchor="middle" font-weight="bold">nokta sayısı N, log</text>')
    o.append(f'<text transform="translate(14,{H/2:.0f}) rotate(-90)" text-anchor="middle" font-weight="bold">süre [s], log</text></svg>')
    open(out, 'w').write('\n'.join(o))

if __name__ == '__main__':
    d, prefix = sys.argv[1], sys.argv[2]
    rows = read(os.path.join(d, 'bench_main.csv'))
    pareto(rows, prefix + '_pareto.svg'); scaling(rows, prefix + '_scaling.svg')
    print('wrote', prefix + '_pareto.svg', prefix + '_scaling.svg')
