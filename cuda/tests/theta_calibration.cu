/**
 * theta calibration: which GPU (LBVH) theta gives the same global L2 force
 * error as a given CPU (octree) theta? Runs entirely on the host (the GPU
 * pipeline is emulated with the functions of bh_core.cuh), so it needs no
 * GPU. Speed comparisons in the benchmark are only meaningful at matched
 * accuracy; rerun this whenever the CPU theta or the distributions change.
 *
 * Usage: bh_theta_calibration [N=8000] [cpu_theta=0.5]
 */
#include "benchmark_framework.cuh"
#include "bh_core.cuh"
#include <vector>
#include <numeric>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <string>
using namespace barnes_hut; using namespace barnes_hut::gpu; using namespace barnes_hut::benchmark;
// host emulation of the GPU pipeline (same order as bh_gpu.cu)
static double gpu_l2(const std::vector<float>& m,const std::vector<float>& x,const std::vector<float>& y,const std::vector<float>& z,
  float theta,float eps2,const std::vector<float>& rx,const std::vector<float>& ry,const std::vector<float>& rz,double& cost){
  int n=x.size(); Box b=empty_box(); PointBox pb{x.data(),y.data(),z.data()}; for(int i=0;i<n;++i) b=merge(b,pb(i));
  std::vector<uint64_t> k(n); for(int i=0;i<n;++i) k[i]=morton_key(x[i],y[i],z[i],b);
  std::vector<int> id(n); std::iota(id.begin(),id.end(),0); std::stable_sort(id.begin(),id.end(),[&](int a,int c){return k[a]<k[c];});
  std::vector<uint64_t> sk(n); std::vector<float> sx(n),sy(n),sz(n),sm(n);
  for(int q=0;q<n;++q){int s=id[q]; sk[q]=k[s]; sx[q]=x[s]; sy[q]=y[s]; sz[q]=z[s]; sm[q]=m[s];}
  std::vector<int> L(n-1),R(n-1),P(2*n-1,-1),V(n-1,0); for(int i=0;i<n-1;++i) build_internal_node(i,sk.data(),n,L.data(),R.data(),P.data());
  std::vector<float4> cm(2*n-1),lo(2*n-1),hi(2*n-1); TreeView t{n,L.data(),R.data(),P.data(),cm.data(),lo.data(),hi.data(),V.data()};
  for(int j=0;j<n;++j) refit_from_leaf(j,sx.data(),sy.data(),sz.data(),sm.data(),t,1.0f/(theta*theta));
  double num=0,den=0; cost=0;
  for(int q=0;q<n;++q){ TraversalCounts c; float3 f=force_on<64>(q,L.data(),R.data(),cm.data(),lo.data(),hi.data(),n,eps2,1.0f,c);
    int o=id[q]; double dx=f.x-rx[o],dy=f.y-ry[o],dz=f.z-rz[o]; num+=dx*dx+dy*dy+dz*dz; den+=double(rx[o])*rx[o]+double(ry[o])*ry[o]+double(rz[o])*rz[o]; cost+=c.particle_particle+c.particle_cell;}
  cost/=n; return std::sqrt(num/den);
}
int main(int argc, char** argv){
  const size_t N = argc > 1 ? std::stoul(argv[1]) : 8000; const float cpu_theta = argc > 2 ? std::stof(argv[2]) : 0.5f;
  BenchmarkConfig cfg;
  std::printf("N=%zu, CPU theta=%.2f, global L2 error vs double direct sum\n", N, cpu_theta);
  for (auto d : {Distribution::Uniform, Distribution::Clustered, Distribution::Plummer, Distribution::Disk}) {
    std::vector<float> m(N),x(N),y(N),z(N),vx(N),vy(N),vz(N),fx(N),fy(N),fz(N),rx(N),ry(N),rz(N);
    generate_distribution(d,N,m.data(),x.data(),y.data(),z.data(),vx.data(),vy.data(),vz.data(),42);
    std::vector<size_t> all(N); std::iota(all.begin(),all.end(),size_t{0});
    reference_forces_direct(m.data(),x.data(),y.data(),z.data(),N,cfg.epsilon_squared,cfg.gravity,all,rx.data(),ry.data(),rz.data());
    for (float th : {cpu_theta}) { cfg.theta=th; cpu_bh_forces(m.data(),x.data(),y.data(),z.data(),N,cfg,fx.data(),fy.data(),fz.data());
      auto r=validate_forces(rx.data(),ry.data(),rz.data(),fx.data(),fy.data(),fz.data(),N,1.0f);
      std::printf("%-10s CPU theta=%.2f L2=%.3e\n",distribution_name(d),th,r.l2_norm_error);}
    for (float th : {0.2f,0.25f,0.3f,0.35f,0.4f,0.5f}) { double cost; double l2=gpu_l2(m,x,y,z,th,cfg.epsilon_squared,rx,ry,rz,cost);
      std::printf("%-10s GPU theta=%.2f L2=%.3e cost=%.0f\n",distribution_name(d),th,l2,cost);}
  }
}
