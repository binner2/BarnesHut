"""Render coloured PLY clouds (from this pipeline) as depth-camera-view images.
Points are projected with the synthetic intrinsics; nearest point wins (z-buffer)."""
import sys, numpy as np
from PIL import Image

def read_ply(path):
    with open(path,'rb') as f:
        props=[]; n=0; binary=True
        while True:
            line=f.readline().decode().strip()
            if line.startswith('format'): binary='binary' in line
            elif line.startswith('element vertex'): n=int(line.split()[2])
            elif line.startswith('property'): 
                t,name=line.split()[1],line.split()[2]; props.append((name,{'float':'<f4','uchar':'u1','double':'<f8','int':'<i4'}[t]))
            elif line=='end_header': break
        dt=np.dtype(props)
        if binary: arr=np.frombuffer(f.read(dt.itemsize*n),dtype=dt)
        else: arr=np.loadtxt(f,dtype=dt,max_rows=n)
    return arr

def render(path, out, fx=525,fy=525,cx=319.5,cy=239.5,w=640,h=480, splat=1, bg=(250,250,248)):
    a=read_ply(path)
    x,y,z=a['x'].astype(float),a['y'].astype(float),a['z'].astype(float)
    if 'red' in a.dtype.names: col=np.stack([a['red'],a['green'],a['blue']],1).astype(np.uint8)
    else: col=np.full((len(x),3),90,np.uint8)
    u=np.round(fx*x/z+cx).astype(int); v=np.round(fy*y/z+cy).astype(int)
    img=np.full((h,w,3),bg,np.uint8); zb=np.full((h,w),np.inf)
    order=np.argsort(-z)  # far first, near overwrite
    for du in range(-splat,splat+1):
        for dv in range(-splat,splat+1):
            uu=u[order]+du; vv=v[order]+dv
            ok=(uu>=0)&(uu<w)&(vv>=0)&(vv<h)
            img[vv[ok],uu[ok]]=col[order][ok]
    Image.fromarray(img).save(out, quality=82, method=6)
    print(out, len(x))

if __name__=='__main__':
    render(*sys.argv[1:3], splat=int(sys.argv[3]) if len(sys.argv)>3 else 1)
