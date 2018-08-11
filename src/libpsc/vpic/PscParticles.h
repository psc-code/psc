
#pragma once

template<class ParticlesBase, class FA, class IA, class AA, class HA>
struct PscParticles : ParticlesBase
{
  typedef ParticlesBase Base;
  typedef PscParticles<Base, FA, IA, AA, HA> Particles;
  typedef FA FieldArray;
  typedef IA Interpolator;
  typedef AA Accumulator;
  typedef HA HydroArray;

  using Base::Base;
  using typename Base::iterator;
  using typename Base::const_iterator;
  using typename Base::Species;
  using typename Base::Particle;
  using typename Base::Grid;

  // ----------------------------------------------------------------------
  // uncenter_p
  
  static void uncenter_p_pipeline(Species* sp, /*const*/ Interpolator& interpolator,
				  int off, int cnt)
  {
    const Grid* g = sp->grid();
    const typename Interpolator::Element* f;
    // For backward half advance
    const float qdt_2mc = -(sp->q * g->dt) / (2*sp->m * g->cvac);
    const float qdt_4mc = 0.5 * qdt_2mc; // For backward half rotate
    const float one = 1.;
    const float one_third = 1./3.;
    const float two_fifteenths = 2./15.;

    float dx, dy, dz, ux, uy, uz;
    float hax, hay, haz, cbx, cby, cbz;
    float v0, v1, v2, v3, v4;

    int ii;

    int n = cnt;
    Particle* p = sp->p + off;

    // Process particles for this pipeline

    for (; n; n--, p++) {
      dx   = p->dx;                            // Load position
      dy   = p->dy;
      dz   = p->dz;
      ii   = p->i;
      f    = &interpolator[ii];                // Interpolate E
      hax  = qdt_2mc*(    ( f->ex    + dy*f->dexdy    ) +
		       dz*( f->dexdz + dy*f->d2exdydz ) );
      hay  = qdt_2mc*(    ( f->ey    + dz*f->deydz    ) +
		       dx*( f->deydx + dz*f->d2eydzdx ) );
      haz  = qdt_2mc*(    ( f->ez    + dx*f->dezdx    ) +
		       dy*( f->dezdy + dx*f->d2ezdxdy ) );
      cbx  = f->cbx + dx*f->dcbxdx;            // Interpolate B
      cby  = f->cby + dy*f->dcbydy;
      cbz  = f->cbz + dz*f->dcbzdz;
      ux   = p->ux;                            // Load momentum
      uy   = p->uy;
      uz   = p->uz;
      v0   = qdt_4mc/(float)sqrt(one + (ux*ux + (uy*uy + uz*uz)));
      /**/                                     // Boris - scalars
      v1   = cbx*cbx + (cby*cby + cbz*cbz);
      v2   = (v0*v0)*v1;
      v3   = v0*(one+v2*(one_third+v2*two_fifteenths));
      v4   = v3/(one+v1*(v3*v3));
      v4  += v4;
      v0   = ux + v3*( uy*cbz - uz*cby );      // Boris - uprime
      v1   = uy + v3*( uz*cbx - ux*cbz );
      v2   = uz + v3*( ux*cby - uy*cbx );
      ux  += v4*( v1*cbz - v2*cby );           // Boris - rotation
      uy  += v4*( v2*cbx - v0*cbz );
      uz  += v4*( v0*cby - v1*cbx );
      ux  += hax;                              // Half advance E
      uy  += hay;
      uz  += haz;
      p->ux = ux;                              // Store momentum
      p->uy = uy;
      p->uz = uz;
    }
  }

#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

  static void uncenter_p_pipeline_v4(Species* sp, /*const*/ Interpolator& interpolator,
				     int off, int cnt)
  {
    using namespace v4;
    const Grid* g = sp->grid();
    const typename Interpolator::Element * ALIGNED(128) f0 = interpolator.data();

    Particle             * ALIGNED(128) p;
    const float          * ALIGNED(16)  vp0;
    const float          * ALIGNED(16)  vp1;
    const float          * ALIGNED(16)  vp2;
    const float          * ALIGNED(16)  vp3;

    const float _qdt_2mc = (sp->q * g->dt) / (2*sp->m * g->cvac);
    
    const v4float qdt_2mc(    -_qdt_2mc); // For backward half advance
    const v4float qdt_4mc(-0.5*_qdt_2mc); // For backward half Boris rotate
    const v4float one(1.);
    const v4float one_third(1./3.);
    const v4float two_fifteenths(2./15.);

    v4float dx, dy, dz, ux, uy, uz, q;
    v4float hax, hay, haz, cbx, cby, cbz;
    v4float v0, v1, v2, v3, v4, v5;
    v4int ii;

    // Determine which particle quads this pipeline processes

    assert(off % 4 == 0);
    assert(cnt % 4 == 0);
    int nq = cnt >> 2;
    p = sp->p + off;

    // Process the particle quads for this pipeline

    for( ; nq; nq--, p+=4 ) {
      load_4x4_tr(&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx,dx,dy,dz,ii);

      // Interpolate fields
      vp0 = (const float * ALIGNED(16))(f0 + ii(0));
      vp1 = (const float * ALIGNED(16))(f0 + ii(1));
      vp2 = (const float * ALIGNED(16))(f0 + ii(2));
      vp3 = (const float * ALIGNED(16))(f0 + ii(3));
      load_4x4_tr(vp0,  vp1,  vp2,  vp3,  hax,v0,v1,v2); hax = qdt_2mc*fma( fma( dy, v2, v1 ), dz, fma( dy, v0, hax ) );
      load_4x4_tr(vp0+4,vp1+4,vp2+4,vp3+4,hay,v3,v4,v5); hay = qdt_2mc*fma( fma( dz, v5, v4 ), dx, fma( dz, v3, hay ) );
      load_4x4_tr(vp0+8,vp1+8,vp2+8,vp3+8,haz,v0,v1,v2); haz = qdt_2mc*fma( fma( dx, v2, v1 ), dy, fma( dx, v0, haz ) );
      load_4x4_tr(vp0+12,vp1+12,vp2+12,vp3+12,cbx,v3,cby,v4); cbx = fma( v3, dx, cbx );
      /**/                                                    cby = fma( v4, dy, cby );
      load_4x2_tr(vp0+16,vp1+16,vp2+16,vp3+16,cbz,v5);        cbz = fma( v5, dz, cbz );

      // Update momentum
      load_4x4_tr(&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux,ux,uy,uz,q);
      /**/                                              // Could use load_4x3_tr
      v0  = qdt_4mc*rsqrt( one + fma( ux,ux, fma( uy,uy, uz*uz ) ) );
      v1  = fma( cbx,cbx, fma( cby,cby, cbz*cbz ) );
      v2  = (v0*v0)*v1;
      v3  = v0*fma( v2, fma( v2, two_fifteenths, one_third ), one );
      v4  = v3*rcp( fma( v3*v3, v1, one ) ); v4 += v4;
      v0  = fma( fms( uy,cbz, uz*cby ), v3, ux );
      v1  = fma( fms( uz,cbx, ux*cbz ), v3, uy );
      v2  = fma( fms( ux,cby, uy*cbx ), v3, uz );
      ux  = fma( fms( v1,cbz, v2*cby ), v4, ux );
      uy  = fma( fms( v2,cbx, v0*cbz ), v4, uy );
      uz  = fma( fms( v0,cby, v1*cbx ), v4, uz );
      ux += hax;
      uy += hay;
      uz += haz;
      store_4x4_tr(ux,uy,uz,q,&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux);
      /**/                                              // Could use store_4x3_tr
    }
  }

#endif
  
  static void uncenter_p(Species* sp, Interpolator& interpolator)
  {
    int cnt = sp->np & ~15;
#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)
    uncenter_p_pipeline_v4(sp, interpolator, 0, cnt);
#else
    uncenter_p_pipeline(sp, interpolator, 0, cnt);
#endif
    uncenter_p_pipeline(sp, interpolator, cnt, sp->np - cnt);
  }

  // ----------------------------------------------------------------------
  // sort_p

  static void sort_p(Species* sp)
  {
    const Grid* g = sp->grid();
    sp->last_sorted = g->step;

    int n_prts = sp->np;
    int vl = VOXEL(1,1,1,             g->nx,g->ny,g->nz);
    int vh = VOXEL(g->nx,g->ny,g->nz, g->nx,g->ny,g->nz) + 1;

    static int * RESTRICT ALIGNED(128) next;
    if (!next) {
      next = new int[g->nv];
    }
    int * RESTRICT ALIGNED(128) partition = sp->partition;

    static Particle * RESTRICT ALIGNED(128) p_aux;
    static size_t n_alloced;
    if (n_prts > n_alloced) {
      delete[] p_aux;
      p_aux = new Particle[n_prts];
      n_alloced = n_prts;
    }
    Particle* RESTRICT ALIGNED(128) p = sp->p;

    // zero counts
    for (int v = vl; v < vh; v++) {
      next[v] = 0;
    }
    
    // find counts
    for (int i = 0; i < n_prts; i++) {
      next[p[i].i]++;
    }
    
    // prefix sum
    int sum = 0;
    for (int v = vl; v < vh; v++) {
      int count = next[v];
      next[v] = sum;
      partition[v] = sum;
      sum += count;
    }
    partition[vh] = sum;
    
    // reorder
    for(int i = 0; i < n_prts; i++) {
      int v = p[i].i;
      int j = next[v]++;
      p_aux[j] = p[i];
    }

    // fix up unused part of partition
    for (int i = 0; i < vl; i++) {
      partition[i] = 0;
    }
    for (int i = vh; i < g->nv; i++) {
      partition[i] = n_prts;
    }

    // OPT: just swap pointer?
    for (int i = 0; i < n_prts; i++) {
      p[i] = p_aux[i];
    }
  }

  // ----------------------------------------------------------------------
  // energy_p

  static double energy_p_pipeline(const_iterator sp,
				  Interpolator &interpolator,
				  int n0, int n1)
  {
    const Grid* g = sp->grid();
    const typename Interpolator::Element * RESTRICT ALIGNED(128) f = interpolator.data();
    const Particle       * RESTRICT ALIGNED(32)  p = sp->p;
    const float qdt_2mc = (sp->q*g->dt)/(2*sp->m*g->cvac);
    const float msp     = sp->m;
    const float one     = 1;

    float dx, dy, dz;
    float v0, v1, v2;

    double en = 0;

    int i, n;

    n1 += n0;

    // Process particles quads for this pipeline

    for( n=n0; n<n1; n++ ) {
      dx  = p[n].dx;
      dy  = p[n].dy;
      dz  = p[n].dz;
      i   = p[n].i;
      v0  = p[n].ux + qdt_2mc*(    ( f[i].ex    + dy*f[i].dexdy    ) +
				   dz*( f[i].dexdz + dy*f[i].d2exdydz ) );
      v1  = p[n].uy + qdt_2mc*(    ( f[i].ey    + dz*f[i].deydz    ) +
				   dx*( f[i].deydx + dz*f[i].d2eydzdx ) );
      v2  = p[n].uz + qdt_2mc*(    ( f[i].ez    + dx*f[i].dezdx    ) +
				   dy*( f[i].dezdy + dx*f[i].d2ezdxdy ) );
      v0  = v0*v0 + v1*v1 + v2*v2;
      v0  = (msp * p[n].w) * (v0 / (one + sqrtf(one + v0)));
      en += (double)v0;
    }

    return en;
  }

#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

  static double energy_p_pipeline_v4(const_iterator sp,
				     Interpolator &interpolator,
				     int n0, int n1)
  {
    using namespace v4;

    const Grid* g = sp->grid();

    const typename Interpolator::Element * RESTRICT ALIGNED(128) f = interpolator.data();

    const float          * RESTRICT ALIGNED(16)  vp0;
    const float          * RESTRICT ALIGNED(16)  vp1;
    const float          * RESTRICT ALIGNED(16)  vp2;
    const float          * RESTRICT ALIGNED(16)  vp3;

    const v4float qdt_2mc((sp->q*g->dt)/(2*sp->m*g->cvac));
    const v4float msp(sp->m);
    const v4float one(1.);

    v4float dx, dy, dz;
    v4float ex, ey, ez;
    v4float v0, v1, v2, w;
    v4int i;

    double en0 = 0, en1 = 0, en2 = 0, en3 = 0;

    // Determine which particle quads this pipeline processes

    const Particle * RESTRICT ALIGNED(32)  p = sp->p + n0;
    int nq = n1 >> 2;

    // Process the particle quads for this pipeline

    for( ; nq; nq--, p+=4 ) {
      load_4x4_tr(&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx,dx,dy,dz,i);

      // Interpolate fields

      vp0 = (float *)(f + i(0));
      vp1 = (float *)(f + i(1));
      vp2 = (float *)(f + i(2));
      vp3 = (float *)(f + i(3));
      load_4x4_tr(vp0,  vp1,  vp2,  vp3,  ex,v0,v1,v2); ex = fma( fma( dy, v2, v1 ), dz, fma( dy, v0, ex ) );
      load_4x4_tr(vp0+4,vp1+4,vp2+4,vp3+4,ey,v0,v1,v2); ey = fma( fma( dz, v2, v1 ), dx, fma( dz, v0, ey ) );
      load_4x4_tr(vp0+8,vp1+8,vp2+8,vp3+8,ez,v0,v1,v2); ez = fma( fma( dx, v2, v1 ), dy, fma( dx, v0, ez ) );

      // Update momentum to half step
      // (note Boris rotation does not change energy so it is unnecessary)

      load_4x4_tr(&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux,v0,v1,v2,w);
      v0  = fma( ex, qdt_2mc, v0 );
      v1  = fma( ey, qdt_2mc, v1 );
      v2  = fma( ez, qdt_2mc, v2 );

      // Accumulate energy

      v0 = fma( v0,v0, fma( v1,v1, v2*v2 ) );
      v0 = (msp * w) * (v0 / (one + sqrt(one + v0))); 
      en0 += (double)v0(0);
      en1 += (double)v0(1);
      en2 += (double)v0(2);
      en3 += (double)v0(3);
    }

    return en0 + en1 + en2 + en3;
  }

#endif

  static double energy_p(const_iterator sp, Interpolator& interpolator)
  {
    const Grid* g = sp->grid();
    int cnt = sp->np & ~15;
#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)
    double local = energy_p_pipeline_v4(sp, interpolator, 0, cnt);
#else
    double local = energy_p_pipeline(sp, interpolator, 0, cnt);
#endif
    local += energy_p_pipeline(sp, interpolator, cnt, sp->np - cnt);

    double global;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, psc_comm_world);
    return global * sqr(g->cvac);
  }

};


