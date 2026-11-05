"""Assemble the report source: inject measured numbers, SVG charts and images
into the template. Run after the benchmark CSVs exist."""
import base64, csv, math, sys, glob, os
TPL = sys.argv[1]; OUT = sys.argv[2]
# results dir (CSV) and image dir; defaults match the repository layout
HERE = os.path.dirname(os.path.abspath(__file__))
S = sys.argv[3] if len(sys.argv) > 3 else os.path.join(HERE, '..', '..', 'experiments', 'results', '2026-11-04')
IMG = sys.argv[4] if len(sys.argv) > 4 else os.path.join(S, 'img')

def read_csv(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k in ('n','theta','kappa','cutoff_sigmas','leaf','threads','build_s','time_s','far_evals','near_evals','max_rel_err','mean_rel_err','rms_rel_err'):
            if k in r and r[k] != '':
                r[k] = float(r[k])
        ex = {}
        for kv in (r.get('extra') or '').split(';'):
            if '=' in kv:
                a,b = kv.split('=',1); ex[a]=float(b)
        r['ex'] = ex
    return rows

main = read_csv(f'{S}/bench_main.csv')
leaf_rows = []
for p in sorted(glob.glob(f'{S}/bench_leaf*.csv')):
    leaf_rows += read_csv(p)

def sel(rows, **kw):
    out = []
    for r in rows:
        ok = True
        for k,v in kw.items():
            if isinstance(v,float) or isinstance(v,int):
                if abs(r[k]-v) > 1e-9: ok=False
            elif r[k] != v: ok=False
        if ok: out.append(r)
    return out

sizes = sorted({int(r['n']) for r in main})
nmax = max(sizes)
V = {}
def fmt(x, d=2):
    s = f'{x:.{d}f}'.replace('.', ',')
    return s
def sci(x):
    if x == 0: return '0'
    e = int(math.floor(math.log10(abs(x)))); m = x/10**e
    return f'{m:.1f}'.replace('.',',') + f'·10<sup>{e}</sup>'
def pct(x, d=1): return fmt(100*x, d) + ' %'

exact = sel(main, algorithm='kde_exact', n=nmax)[0]
V['nmax'] = f'{nmax:,}'.replace(',', ' ')
V['threads'] = str(int(exact['threads']))
V['leaf_main'] = str(int(exact['leaf']))
V['exact_time_nmax'] = fmt(exact['time_s'],1)
V['exact_evals_nmax'] = sci(exact['near_evals'])
for th in (0.3,0.5,0.7,1.0):
    for ka in (0.0,1.0,0.5):
        rs = sel(main, algorithm='kde', n=nmax, theta=th, kappa=ka, cutoff_sigmas=4.0)
        if not rs: continue
        r = rs[0]; tag = f't{int(th*10):02d}_k{int(ka*10):02d}'
        V[f'kde_time_{tag}'] = fmt(r['time_s'],2)
        V[f'kde_rms_{tag}'] = pct(r['rms_rel_err'],2)
        V[f'kde_max_{tag}'] = pct(r['max_rel_err'],1)
        V[f'kde_speedup_{tag}'] = fmt(exact['time_s']/r['time_s'],0)
        V[f'kde_perq_{tag}'] = fmt((r['far_evals']+r['near_evals'])/nmax,0)
# cutoff-only (theta 0, cutoff 4h) at nmax
r = sel(main, algorithm='kde', n=nmax, theta=0.0, cutoff_sigmas=4.0)[0]
V['kde_time_cut'] = fmt(r['time_s'],2); V['kde_rms_cut'] = pct(r['rms_rel_err'],3); V['kde_max_cut']=pct(r['max_rel_err'],2)
V['kde_speedup_cut'] = fmt(exact['time_s']/r['time_s'],0)
V['kde_perq_cut'] = fmt((r['far_evals']+r['near_evals'])/nmax,0)
V['exact_perq'] = fmt(exact['near_evals']/nmax,0)
# theta 0.5 without cutoff
r = sel(main, algorithm='kde', n=nmax, theta=0.5, kappa=0.0, cutoff_sigmas=0.0)[0]
V['kde_time_t05_nocut'] = fmt(r['time_s'],2)
# mean shift
ms_ex = sel(main, algorithm='meanshift', n=nmax, theta=0.0)[0]
V['ms_exact_time'] = fmt(ms_ex['time_s'],1); V['ms_exact_clusters'] = str(int(ms_ex['ex']['clusters']))
for th in (0.3,0.5,0.7,1.0):
    r = sel(main, algorithm='meanshift', n=nmax, theta=th)[0]
    tag = f't{int(th*10):02d}'
    V[f'ms_time_{tag}'] = fmt(r['time_s'],1); V[f'ms_speedup_{tag}'] = fmt(ms_ex['time_s']/r['time_s'],1)
    V[f'ms_agree_{tag}'] = pct(r['ex']['agreement'],1); V[f'ms_clusters_{tag}'] = str(int(r['ex']['clusters']))
# wlop
w_ex = sel(main, algorithm='wlop', n=nmax, theta=0.0)[0]
V['wlop_exact_time'] = fmt(w_ex['time_s'],1); V['wlop_exact_cov'] = fmt(w_ex['ex']['nn_cov'],3); V['wlop_m'] = str(int(w_ex['ex']['m']))
for th in (0.3,0.5,0.7,1.0):
    r = sel(main, algorithm='wlop', n=nmax, theta=th)[0]; tag=f't{int(th*10):02d}'
    V[f'wlop_time_{tag}'] = fmt(r['time_s'],1); V[f'wlop_cov_{tag}'] = fmt(r['ex']['nn_cov'],3); V[f'wlop_speedup_{tag}'] = fmt(w_ex['time_s']/r['time_s'],1)
V['build_time_nmax'] = fmt(exact['build_s']*1000,0)

# ---------- tables ----------
def table_kde():
    rows = []
    for th in (0.3,0.5,0.7,1.0):
        for ka in (0.0,1.0,0.5):
            rs = sel(main, algorithm='kde', n=nmax, theta=th, kappa=ka, cutoff_sigmas=4.0)
            if not rs: continue
            r = rs[0]
            kap = '—' if ka == 0 else fmt(ka,1)
            rows.append(f"<tr><td class='num'>{fmt(th,1)}</td><td class='num'>{kap}</td><td class='num'>{fmt(r['time_s'],2)}</td><td class='num'>{fmt(exact['time_s']/r['time_s'],0)}×</td><td class='num'>{fmt((r['far_evals']+r['near_evals'])/nmax,0)}</td><td class='num'>{pct(r['rms_rel_err'],2)}</td><td class='num'>{pct(r['max_rel_err'],1)}</td></tr>")
    r = sel(main, algorithm='kde', n=nmax, theta=0.0, cutoff_sigmas=4.0)[0]
    rows.insert(0, f"<tr><td class='num'>0 (kesme 4h)</td><td class='num'>—</td><td class='num'>{fmt(r['time_s'],2)}</td><td class='num'>{fmt(exact['time_s']/r['time_s'],0)}×</td><td class='num'>{fmt((r['far_evals']+r['near_evals'])/nmax,0)}</td><td class='num'>{pct(r['rms_rel_err'],3)}</td><td class='num'>{pct(r['max_rel_err'],2)}</td></tr>")
    rows.insert(0, f"<tr><td class='num'>0 (kesin)</td><td class='num'>—</td><td class='num'>{fmt(exact['time_s'],2)}</td><td class='num'>1×</td><td class='num'>{fmt(exact['near_evals']/nmax,0)}</td><td class='num'>0</td><td class='num'>0</td></tr>")
    return '\n'.join(rows)
V['table_kde'] = table_kde()

def table_ms_wlop():
    rows = []
    for th in (0.0,0.3,0.5,0.7,1.0):
        m = sel(main, algorithm='meanshift', n=nmax, theta=th)[0]
        w = sel(main, algorithm='wlop', n=nmax, theta=th)[0]
        agree = '100 %' if th == 0 else pct(m['ex']['agreement'],1)
        label = 'kesin' if th == 0 else fmt(th,1)
        rows.append(f"<tr><td class='num'>{label}</td><td class='num'>{fmt(m['time_s'],1)}</td><td class='num'>{fmt(ms_ex['time_s']/m['time_s'],1)}×</td><td class='num'>{int(m['ex']['clusters'])}</td><td class='num'>{agree}</td><td class='num'>{fmt(w['time_s'],1)}</td><td class='num'>{fmt(w_ex['time_s']/w['time_s'],1)}×</td><td class='num'>{fmt(w['ex']['nn_cov'],3)}</td></tr>")
    return '\n'.join(rows)
V['table_ms_wlop'] = table_ms_wlop()

def table_leaf():
    rows = []
    allrows = leaf_rows + sel(main, n=50000)
    for L in (4,8,16,32):
        ex = sel(allrows, algorithm='kde_exact', leaf=float(L), n=50000)
        if not ex: continue
        ex = ex[0]
        cells = [str(L), fmt(ex['build_s']*1000,0), fmt(ex['time_s'],2)]
        for th,ka in ((0.5,0.0),(0.5,1.0),(0.7,0.0),(0.7,1.0)):
            r = sel(allrows, algorithm='kde', leaf=float(L), n=50000, theta=th, kappa=ka, cutoff_sigmas=4.0)
            if r:
                r = r[0]; cells.append(f"{fmt(r['time_s'],2)} s / {pct(r['rms_rel_err'],2)}")
            else: cells.append('—')
        rows.append('<tr>' + ''.join(f"<td class='num'>{c}</td>" for c in cells) + '</tr>')
    return '\n'.join(rows)
V['table_leaf'] = table_leaf()

# ---------- charts ----------
def svg_pareto():
    W,H = 720,420; ml,mr,mt,mb = 64,20,16,52
    pts = []
    for th in (0.3,0.5,0.7,1.0):
        for ka,cls in ((0.0,'s1'),(1.0,'s2'),(0.5,'s3')):
            rs = sel(main, algorithm='kde', n=nmax, theta=th, kappa=ka, cutoff_sigmas=4.0)
            if rs: pts.append((rs[0]['time_s'], max(rs[0]['rms_rel_err'],1e-6), th, ka, cls))
    cut = sel(main, algorithm='kde', n=nmax, theta=0.0, cutoff_sigmas=4.0)[0]
    pts.append((cut['time_s'], max(cut['rms_rel_err'],1e-6), 0.0, 0.0, 's0'))
    xs = [p[0] for p in pts] + [exact['time_s']]; ys = [p[1] for p in pts]
    x0,x1 = math.log10(min(xs))-0.15, math.log10(max(xs))+0.15
    y0,y1 = math.log10(min(ys))-0.3, math.log10(max(ys))+0.3
    def X(v): return ml + (math.log10(v)-x0)/(x1-x0)*(W-ml-mr)
    def Y(v): return H-mb - (math.log10(v)-y0)/(y1-y0)*(H-mt-mb)
    o = [f'<svg class="chart" role="img" viewBox="0 0 {W} {H}" style="min-width:{int(W*0.92)+4}px"><title>Barnes-Hut KDE: hesap süresi ile göreli RMS hata, n={nmax}</title>']
    # grid
    for e in range(int(math.floor(x0)), int(math.ceil(x1))+1):
        for m in (1,2,5):
            v = m*10**e
            if x0 <= math.log10(v) <= x1:
                o.append(f'<line class="grid" x1="{X(v):.1f}" y1="{mt}" x2="{X(v):.1f}" y2="{H-mb}"/>')
                if m==1 or (x1-x0)<1.5: o.append(f'<text class="tick" x="{X(v):.1f}" y="{H-mb+16}" text-anchor="middle">{fmt(v,2 if v<1 else 0)}</text>')
    for e in range(int(math.floor(y0)), int(math.ceil(y1))+1):
        v = 10**e
        if y0 <= e <= y1:
            o.append(f'<line class="grid" x1="{ml}" y1="{Y(v):.1f}" x2="{W-mr}" y2="{Y(v):.1f}"/>')
            o.append(f'<text class="tick" x="{ml-6}" y="{Y(v)+4:.1f}" text-anchor="end">{pct(v,0) if v>=0.01 else pct(v,2)}</text>')
    o.append(f'<line class="lref" x1="{X(exact["time_s"]):.1f}" y1="{mt}" x2="{X(exact["time_s"]):.1f}" y2="{H-mb}"/>')
    o.append(f'<text class="lab" x="{X(exact["time_s"])-6:.1f}" y="{mt+14}" text-anchor="end">kesin toplam ({fmt(exact["time_s"],1)} s)</text>')
    o.append(f'<text class="lab-b" x="{(ml+W-mr)/2:.0f}" y="{H-10}" text-anchor="middle">hesap süresi [s], log</text>')
    o.append(f'<text class="lab-b" transform="translate(14,{(mt+H-mb)/2:.0f}) rotate(-90)" text-anchor="middle">göreli RMS hata, log</text>')
    for (t,e,th,ka,cls) in pts:
        lbl = ('kesme' if th==0 else f'θ={fmt(th,1)}')
        o.append(f'<circle class="{cls}" cx="{X(t):.1f}" cy="{Y(e):.1f}" r="6"><title>{lbl}{"" if ka==0 else f", κ={fmt(ka,1)}"}: {fmt(t,2)} s, RMS {pct(e,2)}</title></circle>')
        if cls in ('s1','s0'):
            o.append(f'<text class="lab" x="{X(t)+9:.1f}" y="{Y(e)-7:.1f}">{lbl}</text>')
    o.append('</svg>')
    return '\n'.join(o)
V['svg_pareto'] = svg_pareto()

def svg_scaling():
    W,H = 720,400; ml,mr,mt,mb = 64,20,16,52
    series = [('kesin', 's0', lambda n: sel(main, algorithm='kde_exact', n=n)),
              ('θ=0,5', 's1', lambda n: sel(main, algorithm='kde', n=n, theta=0.5, kappa=0.0, cutoff_sigmas=4.0)),
              ('θ=0,5, κ=1', 's2', lambda n: sel(main, algorithm='kde', n=n, theta=0.5, kappa=1.0, cutoff_sigmas=4.0)),
              ('θ=0,7', 's3', lambda n: sel(main, algorithm='kde', n=n, theta=0.7, kappa=0.0, cutoff_sigmas=4.0))]
    data = {name:[(n, f(n)[0]['time_s']) for n in sizes if f(n)] for name,_,f in series}
    allt = [t for v in data.values() for _,t in v]
    x0,x1 = math.log10(min(sizes))-0.1, math.log10(max(sizes))+0.1
    y0,y1 = math.log10(min(allt))-0.2, math.log10(max(allt))+0.2
    def X(v): return ml + (math.log10(v)-x0)/(x1-x0)*(W-ml-mr)
    def Y(v): return H-mb - (math.log10(v)-y0)/(y1-y0)*(H-mt-mb)
    o = [f'<svg class="chart" role="img" viewBox="0 0 {W} {H}" style="min-width:{int(W*0.92)+4}px"><title>KDE hesap süresinin nokta sayısına göre ölçeklenmesi</title>']
    for n in sizes:
        o.append(f'<line class="grid" x1="{X(n):.1f}" y1="{mt}" x2="{X(n):.1f}" y2="{H-mb}"/><text class="tick" x="{X(n):.1f}" y="{H-mb+16}" text-anchor="middle">{n//1000}k</text>')
    for e in range(int(math.floor(y0)), int(math.ceil(y1))+1):
        v=10**e
        if y0<=e<=y1:
            o.append(f'<line class="grid" x1="{ml}" y1="{Y(v):.1f}" x2="{W-mr}" y2="{Y(v):.1f}"/><text class="tick" x="{ml-6}" y="{Y(v)+4:.1f}" text-anchor="end">{fmt(v,2 if v<1 else 0)} s</text>')
    # reference slopes N^2 and N log N anchored at the smallest n of exact
    n0,t0 = data['kesin'][0]
    o.append(f'<line class="lref" x1="{X(n0):.1f}" y1="{Y(t0):.1f}" x2="{X(sizes[-1]):.1f}" y2="{Y(t0*(sizes[-1]/n0)**2):.1f}"/>')
    o.append(f'<text class="lab" x="{X(sizes[-1])-4:.1f}" y="{Y(t0*(sizes[-1]/n0)**2)+14:.1f}" text-anchor="end">∝ N²</text>')
    for name,cls,_ in series:
        v = data[name]
        o.append(f'<polyline class="{"l"+cls[1] if cls!="s0" else "lref"}" fill="none" points="' + ' '.join(f'{X(n):.1f},{Y(t):.1f}' for n,t in v) + '"/>')
        for n,t in v:
            o.append(f'<circle class="{cls}" cx="{X(n):.1f}" cy="{Y(t):.1f}" r="5"><title>{name}, N={n}: {fmt(t,3)} s</title></circle>')
        n,t = v[-1]
        o.append(f'<text class="lab" x="{X(n)+8:.1f}" y="{Y(t)+4:.1f}">{name}</text>')
    o.append(f'<text class="lab-b" x="{(ml+W-mr)/2:.0f}" y="{H-10}" text-anchor="middle">nokta sayısı N, log</text>')
    o.append(f'<text class="lab-b" transform="translate(14,{(mt+H-mb)/2:.0f}) rotate(-90)" text-anchor="middle">süre [s], log</text>')
    o.append('</svg>')
    return '\n'.join(o)
V['svg_scaling'] = svg_scaling()

# images
for name in ('scene_labels','density_full','meanshift_ds','wlop_ds','curvature','ds_1cm'):
    p = f'{IMG}/{name}.webp'
    if os.path.exists(p):
        V[f'img_{name}'] = 'data:image/webp;base64,' + base64.b64encode(open(p,'rb').read()).decode()

tpl = open(TPL).read()
missing = []
import re
def rep(m):
    k = m.group(1)
    if k in V: return V[k]
    missing.append(k); return m.group(0)
out = re.sub(r'\{\{([a-zA-Z0-9_]+)\}\}', rep, tpl)
open(OUT,'w').write(out)
print('missing:', sorted(set(missing)))
print('keys:', len(V))
