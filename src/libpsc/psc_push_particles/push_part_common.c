
#define DIM_X 1
#define DIM_Y 2
#define DIM_Z 4
#define DIM_XY (DIM_X | DIM_Y)
#define DIM_XZ (DIM_X | DIM_Z)
#define DIM_YZ (DIM_Y | DIM_Z)

#define S0X(off) s0x[off+2]
#define S0Y(off) s0y[off+2]
#define S0Z(off) s0z[off+2]
#define S1X(off) s1x[off+2]
#define S1Y(off) s1y[off+2]
#define S1Z(off) s1z[off+2]

static void
do_genc_push_part(int p, struct psc_fields *pf, particle_range_t prts)
{
  creal dt = ppsc->dt;
  creal dqs = .5f * ppsc->coeff.eta * dt;
  creal fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;

  creal dxi = 1.f / ppsc->patch[p].dx[0];
#if (DIM & DIM_X)
  creal s0x[5] = {}, s1x[5];
  creal xl = .5f * dt;
  creal fnqxs = ppsc->patch[p].dx[0] * fnqs / dt;
#endif
#if (DIM & DIM_Y)
  creal s0y[5] = {}, s1y[5];
  creal yl = .5f * dt;
  creal fnqys = ppsc->patch[p].dx[1] * fnqs / dt;
  creal dyi = 1.f / ppsc->patch[p].dx[1];
#endif
#if (DIM & DIM_Z)
  creal s0z[5] = {}, s1z[5];
  creal zl = .5f * dt;
  creal fnqzs = ppsc->patch[p].dx[2] * fnqs / dt;
  creal dzi = 1.f / ppsc->patch[p].dx[2];
#endif

  PARTICLE_ITER_LOOP(prt_iter, prts.begin, prts.end) {
    particle_t *part = particle_iter_deref(prt_iter);

    // x^n, p^n -> x^(n+.5), p^n

    creal root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    creal vxi = part->pxi * root;
    creal vyi = part->pyi * root;
    creal vzi = part->pzi * root;

#if (DIM & DIM_X)
    part->xi += vxi * xl;
#endif
#if (DIM & DIM_Y)
    part->yi += vyi * yl;
#endif
#if (DIM & DIM_Z)
    part->zi += vzi * zl;
#endif

    creal u = part->xi * dxi;
    int lg1 = particle_real_nint(u);
#if (DIM & DIM_X)
    creal h1 = lg1-u;
    creal gmx=.5f*(.5f+h1)*(.5f+h1);
    creal g0x=.75f-h1*h1;
    creal g1x=.5f*(.5f-h1)*(.5f-h1);
#endif
#if (DIM & DIM_Y)
    creal v = part->yi * dyi;
    int lg2 = particle_real_nint(v);
    creal h2 = lg2-v;
    creal gmy=.5f*(.5f+h2)*(.5f+h2);
    creal g0y=.75f-h2*h2;
    creal g1y=.5f*(.5f-h2)*(.5f-h2);
#endif
#if (DIM & DIM_Z)
    creal w = part->zi * dzi;
    int lg3 = particle_real_nint(w);
    creal h3 = lg3-w;
    creal gmz=.5f*(.5f+h3)*(.5f+h3);
    creal g0z=.75f-h3*h3;
    creal g1z=.5f*(.5f-h3)*(.5f-h3);
#endif

#if DIM == DIM_XY
    int lg3 = 0;
#endif
    
    // CHARGE DENSITY FORM FACTOR AT (n+.5)*dt 

#if DIM == DIM_XY
    S0X(-1) = .5f*(1.5f-creal_abs(h1-1.f))*(1.5f-creal_abs(h1-1.f));
    S0X(+0) = .75f-creal_abs(h1)*creal_abs(h1);
    S0X(+1) = .5f*(1.5f-creal_abs(h1+1.f))*(1.5f-creal_abs(h1+1.f));
    S0Y(-1) = .5f*(1.5f-creal_abs(h2-1.f))*(1.5f-creal_abs(h2-1.f));
    S0Y(+0) = .75f-creal_abs(h2)*creal_abs(h2);
    S0Y(+1) = .5f*(1.5f-creal_abs(h2+1.f))*(1.5f-creal_abs(h2+1.f));
#elif DIM == DIM_XZ
    S0X(-1) = gmx;
    S0X(+0) = g0x;
    S0X(+1) = g1x;
    S0Z(-1) = gmz;
    S0Z(+0) = g0z;
    S0Z(+1) = g1z;
#elif DIM == DIM_YZ
    S0Y(-1) = .5f*(1.5f-creal_abs(h2-1.f))*(1.5f-creal_abs(h2-1.f));
    S0Y(+0) = .75f-creal_abs(h2)*creal_abs(h2);
    S0Y(+1) = .5f*(1.5f-creal_abs(h2+1.f))*(1.5f-creal_abs(h2+1.f));
    S0Z(-1) = .5f*(1.5f-creal_abs(h3-1.f))*(1.5f-creal_abs(h3-1.f));
    S0Z(+0) = .75f-creal_abs(h3)*creal_abs(h3);
    S0Z(+1) = .5f*(1.5f-creal_abs(h3+1.f))*(1.5f-creal_abs(h3+1.f));
#endif

#if (DIM & DIM_X)
    u = part->xi * dxi - .5f;
    int lh1 = particle_real_nint(u);
    h1 = lh1 - u;
#else
    u = part->xi * dxi;
    int lh1 = particle_real_nint(u);
#endif
#if (DIM & DIM_Y)
    v = part->yi * dyi - .5f;
    int lh2 = particle_real_nint(v);
    h2 = lh2 - v;
#endif
#if (DIM & DIM_Z)
    w = part->zi * dzi - .5f;
    int lh3 = particle_real_nint(w);
    h3 = lh3 - w;
#else
    int lh3 = 0;
#endif

#if (DIM & DIM_X)
    creal hmx=.5f*(.5f+h1)*(.5f+h1);
    creal h0x=.75f-h1*h1;
    creal h1x=.5f*(.5f-h1)*(.5f-h1);
#endif
#if (DIM & DIM_Y)
    creal hmy=.5f*(.5f+h2)*(.5f+h2);
    creal h0y=.75f-h2*h2;
    creal h1y=.5f*(.5f-h2)*(.5f-h2);
#endif
#if (DIM & DIM_Z)
    creal hmz=.5f*(.5f+h3)*(.5f+h3);
    creal h0z=.75f-h3*h3;
    creal h1z=.5f*(.5f-h3)*(.5f-h3);
#endif

#if DIM == DIM_XY
    // FIELD INTERPOLATION

    creal exq = (gmy*(hmx*F3(pf, EX, lh1-1,lg2-1,lg3) +
		      h0x*F3(pf, EX, lh1  ,lg2-1,lg3) +
		      h1x*F3(pf, EX, lh1+1,lg2-1,lg3)) +
		 g0y*(hmx*F3(pf, EX, lh1-1,lg2  ,lg3) +
		      h0x*F3(pf, EX, lh1  ,lg2  ,lg3) +
		      h1x*F3(pf, EX, lh1+1,lg2  ,lg3)) +
		 g1y*(hmx*F3(pf, EX, lh1-1,lg2+1,lg3) +
		      h0x*F3(pf, EX, lh1  ,lg2+1,lg3) +
		      h1x*F3(pf, EX, lh1+1,lg2+1,lg3)));

    creal eyq = (hmy*(gmx*F3(pf, EY, lg1-1,lh2-1,lg3) +
		      g0x*F3(pf, EY, lg1  ,lh2-1,lg3) +
		      g1x*F3(pf, EY, lg1+1,lh2-1,lg3)) +
		 h0y*(gmx*F3(pf, EY, lg1-1,lh2  ,lg3) +
		      g0x*F3(pf, EY, lg1  ,lh2  ,lg3) +
		      g1x*F3(pf, EY, lg1+1,lh2  ,lg3)) +
		 h1y*(gmx*F3(pf, EY, lg1-1,lh2+1,lg3) +
		      g0x*F3(pf, EY, lg1  ,lh2+1,lg3) +
		      g1x*F3(pf, EY, lg1+1,lh2+1,lg3)));

    creal ezq = (gmy*(gmx*F3(pf, EZ, lg1-1,lg2-1,lh3) +
		      g0x*F3(pf, EZ, lg1  ,lg2-1,lh3) +
		      g1x*F3(pf, EZ, lg1+1,lg2-1,lh3)) +
		 g0y*(gmx*F3(pf, EZ, lg1-1,lg2  ,lh3) +
		      g0x*F3(pf, EZ, lg1  ,lg2  ,lh3) +
		      g1x*F3(pf, EZ, lg1+1,lg2  ,lh3)) +
		 g1y*(gmx*F3(pf, EZ, lg1-1,lg2+1,lh3) +
		      g0x*F3(pf, EZ, lg1  ,lg2+1,lh3) +
		      g1x*F3(pf, EZ, lg1+1,lg2+1,lh3)));

    creal hxq = (hmy*(gmx*F3(pf, HX, lg1-1,lh2-1,lh3) +
		      g0x*F3(pf, HX, lg1  ,lh2-1,lh3) +
		      g1x*F3(pf, HX, lg1+1,lh2-1,lh3)) +
		 h0y*(gmx*F3(pf, HX, lg1-1,lh2  ,lh3) +
		      g0x*F3(pf, HX, lg1  ,lh2  ,lh3) +
		      g1x*F3(pf, HX, lg1+1,lh2  ,lh3)) +
		 h1y*(gmx*F3(pf, HX, lg1-1,lh2+1,lh3) +
		      g0x*F3(pf, HX, lg1  ,lh2+1,lh3) +
		      g1x*F3(pf, HX, lg1+1,lh2+1,lh3)));

    creal hyq = (gmy*(hmx*F3(pf, HY, lh1-1,lg2-1,lh3) +
		      h0x*F3(pf, HY, lh1  ,lg2-1,lh3) +
		      h1x*F3(pf, HY, lh1+1,lg2-1,lh3)) +
		 g0y*(hmx*F3(pf, HY, lh1-1,lg2  ,lh3) +
		      h0x*F3(pf, HY, lh1  ,lg2  ,lh3) +
		      h1x*F3(pf, HY, lh1+1,lg2  ,lh3)) +
		 g1y*(hmx*F3(pf, HY, lh1-1,lg2+1,lh3) +
		      h0x*F3(pf, HY, lh1  ,lg2+1,lh3) +
		      h1x*F3(pf, HY, lh1+1,lg2+1,lh3)));

    creal hzq = (hmy*(hmx*F3(pf, HZ, lh1-1,lh2-1,lg3) +
		      h0x*F3(pf, HZ, lh1  ,lh2-1,lg3) +
		      h1x*F3(pf, HZ, lh1+1,lh2-1,lg3)) +
		 h0y*(hmx*F3(pf, HZ, lh1-1,lh2  ,lg3) +
		      h0x*F3(pf, HZ, lh1  ,lh2  ,lg3) +
		      h1x*F3(pf, HZ, lh1+1,lh2  ,lg3)) +
		 h1y*(hmx*F3(pf, HZ, lh1-1,lh2+1,lg3) +
		      h0x*F3(pf, HZ, lh1  ,lh2+1,lg3) +
		      h1x*F3(pf, HZ, lh1+1,lh2+1,lg3)));
		 
     // c x^(n+.5), p^n -> x^(n+1.0), p^(n+1.0) 

    creal dq = dqs * part->qni / part->mni;
    creal pxm = part->pxi + dq*exq;
    creal pym = part->pyi + dq*eyq;
    creal pzm = part->pzi + dq*ezq;

    root = dq / creal_sqrt(1.f + pxm*pxm + pym*pym + pzm*pzm);
    creal taux = hxq*root;
    creal tauy = hyq*root;
    creal tauz = hzq*root;

    creal tau = 1.f / (1.f + taux*taux + tauy*tauy + tauz*tauz);
    creal pxp = ((1.f+taux*taux-tauy*tauy-tauz*tauz)*pxm + 
		(2.f*taux*tauy+2.f*tauz)*pym + 
		(2.f*taux*tauz-2.f*tauy)*pzm)*tau;
    creal pyp = ((2.f*taux*tauy-2.f*tauz)*pxm +
		(1.f-taux*taux+tauy*tauy-tauz*tauz)*pym +
		(2.f*tauy*tauz+2.f*taux)*pzm)*tau;
    creal pzp = ((2.f*taux*tauz+2.f*tauy)*pxm +
		(2.f*tauy*tauz-2.f*taux)*pym +
		(1.f-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau;
    
    part->pxi = pxp + dq * exq;
    part->pyi = pyp + dq * eyq;
    part->pzi = pzp + dq * ezq;

    root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    vxi = part->pxi * root;
    vyi = part->pyi * root;
    vzi = part->pzi * root;

    part->xi += vxi * xl;
    part->yi += vyi * yl;

    // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)

    creal xi = part->xi + vxi * xl;
    creal yi = part->yi + vyi * yl;

    u = xi * dxi;
    v = yi * dyi;
    int k1 = particle_real_nint(u);
    int k2 = particle_real_nint(v);
    h1 = k1 - u;
    h2 = k2 - v;

    for (int i = -2; i <= 2; i++) {
      S1X(i) = 0.f;
      S1Y(i) = 0.f;
    }

    S1X(k1-lg1-1) = .5f*(1.5f-creal_abs(h1-1.f))*(1.5f-creal_abs(h1-1.f));
    S1X(k1-lg1+0) = .75f-creal_abs(h1)*creal_abs(h1);
    S1X(k1-lg1+1) = .5f*(1.5f-creal_abs(h1+1.f))*(1.5f-creal_abs(h1+1.f));
    S1Y(k2-lg2-1) = .5f*(1.5f-creal_abs(h2-1.f))*(1.5f-creal_abs(h2-1.f));
    S1Y(k2-lg2+0) = .75f-creal_abs(h2)*creal_abs(h2);
    S1Y(k2-lg2+1) = .5f*(1.5f-creal_abs(h2+1.f))*(1.5f-creal_abs(h2+1.f));

    // CURRENT DENSITY AT (n+1.0)*dt

    for (int i = -1; i <= 1; i++) {
      S1X(i) -= S0X(i);
      S1Y(i) -= S0Y(i);
    }

    int l1min, l2min, l1max, l2max;
    
    if (k1 == lg1) {
      l1min = -1; l1max = +1;
    } else if (k1 == lg1 - 1) {
      l1min = -2; l1max = +1;
    } else { // (k1 == lg1 + 1)
      l1min = -1; l1max = +2;
    }

    if (k2 == lg2) {
      l2min = -1; l2max = +1;
    } else if (k2 == lg2 - 1) {
      l2min = -2; l2max = +1;
    } else { // (k2 == lg2 + 1)
      l2min = -1; l2max = +2;
    }

    creal fnqx = part->qni * part->wni * fnqxs;
    creal fnqy = part->qni * part->wni * fnqys;
    creal fnqz = vzi * part->qni * part->wni * fnqs;
    for (int l2 = l2min; l2 <= l2max; l2++) {
      creal jxh = 0.f;
      for (int l1 = l1min; l1 <= l1max; l1++) {
	creal wx = S1X(l1) * (S0Y(l2) + .5f*S1Y(l2));
	creal wz = S0X(l1) * S0Y(l2)
	  + .5f * S1X(l1) * S0Y(l2)
	  + .5f * S0X(l1) * S1Y(l2)
	  + (1.f/3.f) * S1X(l1) * S1Y(l2);

	jxh -= fnqx*wx;
	F3(pf, JXI, lg1+l1,lg2+l2,lg3) += jxh;
	F3(pf, JZI, lg1+l1,lg2+l2,lg3) += fnqz * wz;
      }
    }
    for (int l1 = l1min; l1 <= l1max; l1++) {
      creal jyh = 0.f;
      for (int l2 = l2min; l2 <= l2max; l2++) {
	creal wy = S1Y(l2) * (S0X(l1) + .5f*S1X(l1));

	jyh -= fnqy*wy;
	F3(pf, JYI, lg1+l1,lg2+l2,lg3) += jyh;
      }
    }

#elif DIM == DIM_XZ

    // FIELD INTERPOLATION

#define IP_2ND(pf, m, gx, gy, gz)					\
    (gz##mz*(gx##mx*F3(pf, m, l##gx##1-1,0,l##gz##3-1) +			\
	     gx##0x*F3(pf, m, l##gx##1  ,0,l##gz##3-1) +		\
	     gx##1x*F3(pf, m, l##gx##1+1,0,l##gz##3-1)) +			\
     gz##0z*(gx##mx*F3(pf, m, l##gx##1-1,0,l##gz##3  ) +			\
	     gx##0x*F3(pf, m, l##gx##1  ,0,l##gz##3  ) +			\
	     gx##1x*F3(pf, m, l##gx##1+1,0,l##gz##3  )) +			\
     gz##1z*(gx##mx*F3(pf, m, l##gx##1-1,0,l##gz##3+1) +			\
	     gx##0x*F3(pf, m, l##gx##1  ,0,l##gz##3+1) +			\
	     gx##1x*F3(pf, m, l##gx##1+1,0,l##gz##3+1)))			\
      
    creal exq = IP_2ND(pf, EX, h, g, g);
    creal eyq = IP_2ND(pf, EY, g, h, g);
    creal ezq = IP_2ND(pf, EZ, g, g, h);

    creal hxq = IP_2ND(pf, HX, g, h, h);
    creal hyq = IP_2ND(pf, HY, h, g, h);
    creal hzq = IP_2ND(pf, HZ, h, h, g);

     // c x^(n+.5), p^n -> x^(n+1.0), p^(n+1.0) 

    creal dq = dqs * part->qni / part->mni;
    creal pxm = part->pxi + dq*exq;
    creal pym = part->pyi + dq*eyq;
    creal pzm = part->pzi + dq*ezq;

    root = dq / creal_sqrt(1.f + pxm*pxm + pym*pym + pzm*pzm);
    creal taux = hxq*root;
    creal tauy = hyq*root;
    creal tauz = hzq*root;

    creal tau = 1.f / (1.f + taux*taux + tauy*tauy + tauz*tauz);
    creal pxp = ((1.f+taux*taux-tauy*tauy-tauz*tauz)*pxm + 
		(2.f*taux*tauy+2.f*tauz)*pym + 
		(2.f*taux*tauz-2.f*tauy)*pzm)*tau;
    creal pyp = ((2.f*taux*tauy-2.f*tauz)*pxm +
		(1.f-taux*taux+tauy*tauy-tauz*tauz)*pym +
		(2.f*tauy*tauz+2.f*taux)*pzm)*tau;
    creal pzp = ((2.f*taux*tauz+2.f*tauy)*pxm +
		(2.f*tauy*tauz-2.f*taux)*pym +
		(1.f-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau;
    
    part->pxi = pxp + dq * exq;
    part->pyi = pyp + dq * eyq;
    part->pzi = pzp + dq * ezq;

    root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    vxi = part->pxi * root;
    vyi = part->pyi * root;
    vzi = part->pzi * root;

    part->xi += vxi * xl;
    part->zi += vzi * zl;

    // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)

    creal xi = part->xi + vxi * xl;
    creal zi = part->zi + vzi * zl;

    u = xi * dxi;
    w = zi * dzi;
    int k1 = particle_real_nint(u);
    int k3 = particle_real_nint(w);
    h1 = k1 - u;
    h3 = k3 - w;

    for (int i = -2; i <= 2; i++) {
      S1X(i) = 0.f;
      S1Z(i) = 0.f;
    }

    S1X(k1-lg1-1) = .5f*(.5f+h1)*(.5f+h1);
    S1X(k1-lg1+0) = .75f-h1*h1;
    S1X(k1-lg1+1) = .5f*(.5f-h1)*(.5f-h1);
    S1Z(k3-lg3-1) = .5f*(.5f+h3)*(.5f+h3);
    S1Z(k3-lg3+0) = .75f-h3*h3;
    S1Z(k3-lg3+1) = .5f*(.5f-h3)*(.5f-h3);

    // CURRENT DENSITY AT (n+1.0)*dt

    for (int i = -1; i <= 1; i++) {
      S1X(i) -= S0X(i);
      S1Z(i) -= S0Z(i);
    }

    int l1min, l3min, l1max, l3max;
    
    if (k1 == lg1) {
      l1min = -1; l1max = +1;
    } else if (k1 == lg1 - 1) {
      l1min = -2; l1max = +1;
    } else { // (k1 == lg1 + 1)
      l1min = -1; l1max = +2;
    }

    if (k3 == lg3) {
      l3min = -1; l3max = +1;
    } else if (k3 == lg3 - 1) {
      l3min = -2; l3max = +1;
    } else { // (k3 == lg3 + 1)
      l3min = -1; l3max = +2;
    }

    creal fnqx = part->qni * part->wni * fnqxs;
    for (int l3 = l3min; l3 <= l3max; l3++) {
      creal jxh = 0.f;
      for (int l1 = l1min; l1 < l1max; l1++) {
	creal wx = S1X(l1) * (S0Z(l3) + .5f*S1Z(l3));
	jxh -= fnqx*wx;
	F3(pf, JXI, lg1+l1,0,lg3+l3) += jxh;
      }
    }

    creal fnqy = vyi * part->qni * part->wni * fnqs;
    for (int l3 = l3min; l3 <= l3max; l3++) {
      for (int l1 = l1min; l1 <= l1max; l1++) {
	creal wy = S0X(l1) * S0Z(l3)
	  + .5f * S1X(l1) * S0Z(l3)
	  + .5f * S0X(l1) * S1Z(l3)
	  + (1.f/3.f) * S1X(l1) * S1Z(l3);
	creal jyh = fnqy*wy;
	F3(pf, JYI, lg1+l1,0,lg3+l3) += jyh;
      }
    }

    creal fnqz = part->qni * part->wni * fnqzs;
    for (int l1 = l1min; l1 <= l1max; l1++) {
      creal jzh = 0.f;
      for (int l3 = l3min; l3 < l3max; l3++) {
	creal wz = S1Z(l3) * (S0X(l1) + .5f*S1X(l1));
	jzh -= fnqz*wz;
	F3(pf, JZI, lg1+l1,0,lg3+l3) += jzh;
      }
    }

#elif DIM == DIM_YZ

    // FIELD INTERPOLATION

#define IP_2ND(pf, m, gx, gy, gz)					\
    (gz##mz*(gy##my*F3(pf, m, l##gx##1,l##gy##2-1,l##gz##3-1) +		\
	     gy##0y*F3(pf, m, l##gx##1,l##gy##2  ,l##gz##3-1) +		\
	     gy##1y*F3(pf, m, l##gx##1,l##gy##2+1,l##gz##3-1)) +	\
     gz##0z*(gy##my*F3(pf, m, l##gx##1,l##gy##2-1,l##gz##3  ) +		\
	     gy##0y*F3(pf, m, l##gx##1,l##gy##2  ,l##gz##3  ) +		\
	     gy##1y*F3(pf, m, l##gx##1,l##gy##2+1,l##gz##3  )) +	\
     gz##1z*(gy##my*F3(pf, m, l##gx##1,l##gy##2-1,l##gz##3+1) +		\
	     gy##0y*F3(pf, m, l##gx##1,l##gy##2  ,l##gz##3+1) +		\
	     gy##1y*F3(pf, m, l##gx##1,l##gy##2+1,l##gz##3+1)))		\
      
    creal exq = IP_2ND(pf, EX, h, g, g);
    creal eyq = IP_2ND(pf, EY, g, h, g);
    creal ezq = IP_2ND(pf, EZ, g, g, h);
    creal hxq = IP_2ND(pf, HX, g, h, h);
    creal hyq = IP_2ND(pf, HY, h, g, h);
    creal hzq = IP_2ND(pf, HZ, h, h, g);

     // c x^(n+.5), p^n -> x^(n+1.0), p^(n+1.0) 

    creal dq = dqs * part->qni / part->mni;
    creal pxm = part->pxi + dq*exq;
    creal pym = part->pyi + dq*eyq;
    creal pzm = part->pzi + dq*ezq;

    root = dq / creal_sqrt(1.f + pxm*pxm + pym*pym + pzm*pzm);
    creal taux = hxq*root;
    creal tauy = hyq*root;
    creal tauz = hzq*root;

    creal tau = 1.f / (1.f + taux*taux + tauy*tauy + tauz*tauz);
    creal pxp = ((1.f+taux*taux-tauy*tauy-tauz*tauz)*pxm + 
		(2.f*taux*tauy+2.f*tauz)*pym + 
		(2.f*taux*tauz-2.f*tauy)*pzm)*tau;
    creal pyp = ((2.f*taux*tauy-2.f*tauz)*pxm +
		(1.f-taux*taux+tauy*tauy-tauz*tauz)*pym +
		(2.f*tauy*tauz+2.f*taux)*pzm)*tau;
    creal pzp = ((2.f*taux*tauz+2.f*tauy)*pxm +
		(2.f*tauy*tauz-2.f*taux)*pym +
		(1.f-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau;
    
    part->pxi = pxp + dq * exq;
    part->pyi = pyp + dq * eyq;
    part->pzi = pzp + dq * ezq;

    root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    vxi = part->pxi * root;
    vyi = part->pyi * root;
    vzi = part->pzi * root;

    part->yi += vyi * yl;
    part->zi += vzi * zl;

    // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)

    creal yi = part->yi + vyi * yl;
    creal zi = part->zi + vzi * zl;

    v = yi * dyi;
    w = zi * dzi;
    int k2 = particle_real_nint(v);
    int k3 = particle_real_nint(w);
    h2 = k2 - v;
    h3 = k3 - w;

    for (int i = -2; i <= 2; i++) {
      S1Y(i) = 0.f;
      S1Z(i) = 0.f;
    }

    S1Y(k2-lg2-1) = .5f*(1.5f-creal_abs(h2-1.f))*(1.5f-creal_abs(h2-1.f));
    S1Y(k2-lg2+0) = .75f-creal_abs(h2)*creal_abs(h2);
    S1Y(k2-lg2+1) = .5f*(1.5f-creal_abs(h2+1.f))*(1.5f-creal_abs(h2+1.f));
    S1Z(k3-lg3-1) = .5f*(1.5f-creal_abs(h3-1.f))*(1.5f-creal_abs(h3-1.f));
    S1Z(k3-lg3+0) = .75f-creal_abs(h3)*creal_abs(h3);
    S1Z(k3-lg3+1) = .5f*(1.5f-creal_abs(h3+1.f))*(1.5f-creal_abs(h3+1.f));

    // CURRENT DENSITY AT (n+1.0)*dt

    for (int i = -1; i <= 1; i++) {
      S1Y(i) -= S0Y(i);
      S1Z(i) -= S0Z(i);
    }

    int l2min, l3min, l2max, l3max;
    
    if (k2 == lg2) {
      l2min = -1; l2max = +1;
    } else if (k2 == lg2 - 1) {
      l2min = -2; l2max = +1;
    } else { // (k2 == lg2 + 1)
      l2min = -1; l2max = +2;
    }

    if (k3 == lg3) {
      l3min = -1; l3max = +1;
    } else if (k3 == lg3 - 1) {
      l3min = -2; l3max = +1;
    } else { // (k3 == lg3 + 1)
      l3min = -1; l3max = +2;
    }

    creal jxh;
    creal jyh;
    creal jzh[5];

#define JZH(i) jzh[i+2]

    creal fnqx = vxi * part->qni * part->wni * fnqs;
    creal fnqy = part->qni * part->wni * fnqys;
    creal fnqz = part->qni * part->wni * fnqzs;
    for (int l2 = l2min; l2 <= l2max; l2++) {
      JZH(l2) = 0.f;
    }
    for (int l3 = l3min; l3 <= l3max; l3++) {
      jyh = 0.f;
      for (int l2 = l2min; l2 <= l2max; l2++) {
	creal wx = S0Y(l2) * S0Z(l3)
	  + .5f * S1Y(l2) * S0Z(l3)
	  + .5f * S0Y(l2) * S1Z(l3)
	  + (1.f/3.f) * S1Y(l2) * S1Z(l3);
	creal wy = S1Y(l2) * (S0Z(l3) + .5f*S1Z(l3));
	creal wz = S1Z(l3) * (S0Y(l2) + .5f*S1Y(l2));

	jxh = fnqx*wx;
	jyh -= fnqy*wy;
	JZH(l2) -= fnqz*wz;

	F3(pf, JXI, lg1,lg2+l2,lg3+l3) += jxh;
	F3(pf, JYI, lg1,lg2+l2,lg3+l3) += jyh;
	F3(pf, JZI, lg1,lg2+l2,lg3+l3) += JZH(l2);
      }
    }
#endif
  }
}


