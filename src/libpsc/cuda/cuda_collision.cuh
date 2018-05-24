
#pragma once

#include "cuda_mparticles.h"
#include "cuda_mparticles_sort.cuh"

#include <curand_kernel.h>

#define THREADS_PER_BLOCK 128

// ----------------------------------------------------------------------
// k_curand_setup

__global__
static void k_curand_setup(curandState *d_curand_states)
{
  int bid = blockIdx.x;
  int id = threadIdx.x + bid * THREADS_PER_BLOCK;

  curand_init(1234, id % 1024, 0, &d_curand_states[id]); // FIXME, % 1024 hack
}


// FIXME, replicated from CPU
__device__
static float bc(DMparticlesCuda<BS144> dmprts, float nudt1, int n1, int n2,
		curandState& curand_state)
{
  using real_t = float;
  
  real_t nudt;
  
  real_t pn1,pn2,pn3,pn4;
  real_t p01,p02,pc01,pc02,pc03,pc04;//p03,p04
  real_t px1,py1,pz1,pcx1,pcy1,pcz1;
  real_t px2,py2,pz2,pcx2,pcy2,pcz2;
  real_t px3,py3,pz3,pcx3,pcy3,pcz3;
  real_t px4,py4,pz4,pcx4,pcy4,pcz4;
  real_t h1,h2,h3,h4,ppc,qqc,ss;
  
  real_t m1,m2,m3,m4;
  real_t q1,q2;//,q3,q4;
  //  real_t w1,w2,w3,w4,ww;
  real_t vcx,vcy,vcz;
  real_t bet,gam;
  
  real_t psi,nu;
  real_t nx,ny,nz,nnorm;
  real_t nn1,nx1,ny1,nz1;
  real_t nn2,nx2,ny2,nz2;
  real_t nn3,nx3,ny3,nz3;
  real_t vcx1,vcy1,vcz1;
  real_t vcx2,vcy2,vcz2;
  real_t vcn,vcxr,vcyr,vczr,vcr;
  real_t m12,q12;
  real_t ran1,ran2;

  d_particle prt1, prt2;
  LOAD_PARTICLE_POS(prt1, dmprts.xi4_ , n1);
  LOAD_PARTICLE_POS(prt2, dmprts.xi4_ , n2);
  int kind1 = __float_as_int(prt1.kind_as_float);
  int kind2 = __float_as_int(prt2.kind_as_float);
  q1 = dmprts.q(kind1);
  q2 = dmprts.q(kind2);
  m1 = dmprts.m(kind1);
  m2 = dmprts.m(kind2);

  if (q1*q2 == 0.) {
    return 0.; // no Coulomb collisions with neutrals
  }

  LOAD_PARTICLE_MOM(prt1, dmprts.pxi4_, n1);
  LOAD_PARTICLE_MOM(prt2, dmprts.pxi4_, n2);

  px1 = prt1.pxi[0];
  py1 = prt1.pxi[1];
  pz1 = prt1.pxi[2];

  px2 = prt2.pxi[0];
  py2 = prt2.pxi[1];
  pz2 = prt2.pxi[2];

  px1=m1*px1;
  py1=m1*py1;
  pz1=m1*pz1;
  px2=m2*px2;
  py2=m2*py2;
  pz2=m2*pz2;
  
  
  // determine absolute value of pre-collision momentum in cm-frame
    
  p01=std::sqrt(m1*m1+px1*px1+py1*py1+pz1*pz1);
  p02=std::sqrt(m2*m2+px2*px2+py2*py2+pz2*pz2);
  h1=p01*p02-px1*px2-py1*py2-pz1*pz2;
  ss=m1*m1+m2*m2+2.0*h1;
  h2=ss-m1*m1-m2*m2;
  h3=(h2*h2-4.0*m1*m1*m2*m2)/(4.0*ss);
  if (h3 < 0.) {
    printf("WARNING: ss %g (m1+m2)^2 %g h3 %g in %s\n", // FIXME
	   ss, (m1+m2)*(m1+m2), h3, __FILE__);
    return 0.; // nudt = 0 because no collision
  }
  ppc=std::sqrt(h3);
  
  
  // determine cm-velocity
    
  vcx=(px1+px2)/(p01+p02);
  vcy=(py1+py2)/(p01+p02);
  vcz=(pz1+pz2)/(p01+p02);
  
  nnorm=std::sqrt(vcx*vcx+vcy*vcy+vcz*vcz);
  if (nnorm>0.0) {
    nx=vcx/nnorm;
    ny=vcy/nnorm;
    nz=vcz/nnorm;
  } else {
    nx=0.0;
    ny=0.0;
    nz=0.0;
  }
  bet=nnorm;
  gam=1.0/std::sqrt(1.0-bet*bet);
  
  
  // determine pre-collision momenta in cm-frame
  
      
  pn1=px1*nx+py1*ny+pz1*nz;
  pn2=px2*nx+py2*ny+pz2*nz;
  pc01=std::sqrt(m1*m1+ppc*ppc);
  pcx1=px1+(gam-1.0)*pn1*nx-gam*vcx*p01;
  pcy1=py1+(gam-1.0)*pn1*ny-gam*vcy*p01;
  pcz1=pz1+(gam-1.0)*pn1*nz-gam*vcz*p01;
  pc02=std::sqrt(m2*m2+ppc*ppc);
  pcx2=px2+(gam-1.0)*pn2*nx-gam*vcx*p02;
  pcy2=py2+(gam-1.0)*pn2*ny-gam*vcy*p02;
  pcz2=pz2+(gam-1.0)*pn2*nz-gam*vcz*p02;
      
      
  //  introduce right-handed coordinate system
      
      
  nn1=std::sqrt(pcx1*pcx1+pcy1*pcy1+pcz1*pcz1);
  nn2=std::sqrt(pcx1*pcx1+pcy1*pcy1);
  nn3=nn1*nn2;
  
  if (nn2 != 0.0) {
    nx1=pcx1/nn1;
    ny1=pcy1/nn1;
    nz1=pcz1/nn1;
    
    nx2=pcy1/nn2;
    ny2=-pcx1/nn2;
    nz2=0.0;
    
    nx3=-pcx1*pcz1/nn3;
    ny3=-pcy1*pcz1/nn3;
    nz3=nn2*nn2/nn3;
  } else {
    nx1=0.0;
    ny1=0.0;
    nz1=1.0;
    
    nx2=0.0;
    ny2=1.0;
    nz2=0.0;
    
    nx3=1.0;
    ny3=0.0;
    nz3=0.0;
  }
    
          
  // determine relative particle velocity in cm-frame
          
          
  vcx1=pcx1/pc01;
  vcy1=pcy1/pc01;
  vcz1=pcz1/pc01;
  vcx2=pcx2/pc02;
  vcy2=pcy2/pc02;
  vcz2=pcz2/pc02;
  
  vcn=1.0/(1.0-(vcx1*vcx2+vcy1*vcy2+vcz1*vcz2));
  vcxr=vcn*(vcx1-vcx2);
  vcyr=vcn*(vcy1-vcy2);
  vczr=vcn*(vcz1-vcz2);
  vcr = std::sqrt(vcxr*vcxr+vcyr*vcyr+vczr*vczr);
  if (vcr < 1.0e-20) {
    vcr=1.0e-20;
  }
  
  m3=m1;
  m4=m2;
  /* q3=q1; */
  /* q4=q2; */
  //  c      w3=w1
  //  c      w4=w2
          
          
  // determine absolute value of post-collision momentum in cm-frame
          
  h2=ss-m3*m3-m4*m4;
  h3=(h2*h2-4.0*m3*m3*m4*m4)/(4.0*ss);
  if (h3 < 0.) {
    printf("WARNING: ss %g (m3+m4)^2 %g in %s\n",
	   ss, (m3+m4)*(m3+m4), __FILE__); // FIXME
    return 0.; // nudt = 0 because no collision
  }
          
  qqc=std::sqrt(h3);
  m12=m1*m2/(m1+m2);
  q12=q1*q2;
    
  nudt=nudt1 * q12*q12/(m12*m12*vcr*vcr*vcr);
  
  // event generator of angles for post collision vectors
  
  ran1 = curand_uniform(&curand_state);
  ran2 = curand_uniform(&curand_state);
  if (ran2 < 1e-20) {
    ran2 = 1e-20;
  }
  
  nu = 2.f * M_PI * ran1;
  
  if(nudt<1.0) {                   // small angle collision
    psi=2.f*std::atan(std::sqrt(-.5f*nudt*std::log(ran2)));
  } else {
    psi=std::acos(1.f-2.f*ran2);          // isotropic angles
  }
            
  // determine post-collision momentum in cm-frame
  h1=std::cos(psi);
  h2=std::sin(psi);
  h3=std::sin(nu);
  h4=std::cos(nu);
  
  pc03=std::sqrt(m3*m3+qqc*qqc);
  pcx3=qqc*(h1*nx1+h2*h3*nx2+h2*h4*nx3);
  pcy3=qqc*(h1*ny1+h2*h3*ny2+h2*h4*ny3);
  pcz3=qqc*(h1*nz1+h2*h3*nz2+h2*h4*nz3);
  
  pc04=std::sqrt(m4*m4+qqc*qqc);
  //  c      pcx4=-qqc*(h1*nx1+h2*h3*nx2+h2*h4*nx3)
  //  c      pcy4=-qqc*(h1*ny1+h2*h3*ny2+h2*h4*ny3)
  //  c      pcz4=-qqc*(h1*nz1+h2*h3*nz2+h2*h4*nz3)
  pcx4=-pcx3;
  pcy4=-pcy3;
  pcz4=-pcz3;
  
  // determine post-collision momentum in lab-frame
  
  pn3=pcx3*nx+pcy3*ny+pcz3*nz;
  pn4=pcx4*nx+pcy4*ny+pcz4*nz;
  /* p03=gam*(pc03+bet*pn3); */
  px3=pcx3+(gam-1.0)*pn3*nx+gam*vcx*pc03;
  py3=pcy3+(gam-1.0)*pn3*ny+gam*vcy*pc03;
  pz3=pcz3+(gam-1.0)*pn3*nz+gam*vcz*pc03;
  /* p04=gam*(pc04+bet*pn4); */
  px4=pcx4+(gam-1.0)*pn4*nx+gam*vcx*pc04;
  py4=pcy4+(gam-1.0)*pn4*ny+gam*vcy*pc04;
  pz4=pcz4+(gam-1.0)*pn4*nz+gam*vcz*pc04;
  
  px3=px3/m3;
  py3=py3/m3;
  pz3=pz3/m3;
  px4=px4/m4;
  py4=py4/m4;
  pz4=pz4/m4;
  
  prt1.pxi[0] = px3;
  prt1.pxi[1] = py3;
  prt1.pxi[2] = pz3;
  prt2.pxi[0] = px4;
  prt2.pxi[1] = py4;
  prt2.pxi[2] = pz4;

  STORE_PARTICLE_MOM(prt1, dmprts.pxi4_, n1);
  STORE_PARTICLE_MOM(prt2, dmprts.pxi4_, n2);
  
  return nudt;
}

__global__ static void
k_collide(DMparticlesCuda<BS144> dmprts, uint* d_off, uint* d_id, float nudt0,
	  curandState* d_curand_states, uint n_cells)
{
  using real_t = float;

  int id = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  /* Copy state to local memory for efficiency */
  curandState local_state = d_curand_states[id];

  for (uint bidx = blockIdx.x; bidx < n_cells; bidx += gridDim.x) {
    uint beg = d_off[bidx];
    uint end = d_off[bidx + 1];
    real_t nudt1 = nudt0 * (end - beg & ~1); // somewhat counteract that we don't collide the last particle if odd
    for (uint n = beg + 2*threadIdx.x; n + 1 < end; n += 2*THREADS_PER_BLOCK) {
      //printf("%d/%d: n = %d off %d\n", blockIdx.x, threadIdx.x, n, d_off[blockIdx.x]);
      bc(dmprts, nudt1, d_id[n], d_id[n + 1], local_state);
    }
  }
  
  d_curand_states[id] = local_state;
}

// ======================================================================
// cuda_collision

template<typename cuda_mparticles>
struct cuda_collision
{
  using real_t = typename cuda_mparticles::real_t;
  
  cuda_collision(int interval, double nu, int nicell, double dt)
    : interval_{interval}, nu_{nu}, nicell_(nicell), dt_(dt)
  {}
  
  void operator()(cuda_mparticles& cmprts)
  {
    auto sort_by_cell = cuda_mparticles_sort{cmprts.n_cells()};
    
    sort_by_cell.find_indices_ids(cmprts);
    sort_by_cell.stable_sort_cidx();
    sort_by_cell.find_offsets();
    // for (int c = 0; c <= cmprts.n_cells(); c++) {
    //   printf("off[%d] = %d\n", c, int(sort_by_cell.d_off[c]));
    // }

    static bool first_time = true;
    if (first_time) {
      int blocks = cmprts.n_cells();
      if (blocks > 32768) blocks = 32768;
      dim3 dimGrid(cmprts.n_cells());
      int n_threads = dimGrid.x * THREADS_PER_BLOCK;
      
      cudaError_t ierr;
      ierr = cudaMalloc(&d_curand_states_, n_threads * sizeof(*d_curand_states_));
      cudaCheck(ierr);
      
      k_curand_setup<<<dimGrid, THREADS_PER_BLOCK>>>(d_curand_states_);
      cuda_sync_if_enabled();
      
      first_time = false;
    }
    
    // all particles need to have same weight!
    real_t wni = 1.; // FIXME, there should at least be some assert to enforce this //prts.prt_wni(prts[n_start]);
    real_t nudt0 = wni / nicell_ * interval_ * dt_ * nu_;

    dim3 dimGrid(cmprts.n_cells());
    k_collide<<<dimGrid, THREADS_PER_BLOCK>>>(cmprts, sort_by_cell.d_off.data().get(),
					      sort_by_cell.d_id.data().get(),
					      nudt0, d_curand_states_,
					      cmprts.n_cells());
    cuda_sync_if_enabled();
  }

private:
  int interval_;
  double nu_;
  int nicell_;
  double dt_;
  curandState* d_curand_states_;
};

