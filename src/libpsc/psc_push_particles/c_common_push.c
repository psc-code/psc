
static inline void
calc_vxi(creal vxi[3], particle_t *part)
{
  creal root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
  vxi[0] = part->pxi * root;
  vxi[1] = part->pyi * root;
  vxi[2] = part->pzi * root;
}

static inline void
push_xi(particle_t *part, creal vxi[3], creal dt)
{
  part->yi += vxi[1] * dt;
  part->zi += vxi[2] * dt;
}

static inline void
push_pxi(particle_t *part, creal exq, creal eyq, creal ezq,
	 creal hxq, creal hyq, creal hzq, creal dqs)
{
  creal dq = dqs * part->qni / part->mni;
  creal pxm = part->pxi + dq*exq;
  creal pym = part->pyi + dq*eyq;
  creal pzm = part->pzi + dq*ezq;
  
  creal root = dq / creal_sqrt(1.f + pxm*pxm + pym*pym + pzm*pzm);
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
}

static inline void
find_idx_off_1st(creal xi[3], int lg[3], creal og[3], creal shift,
		 double xb[3], creal dxi[3])
{
  for (int d = 0; d < 3; d++) {
    creal pos = (xi[d] - xb[d]) * dxi[d] + shift;
    lg[d] = fint(pos);
    og[d] = pos - lg[d];
  }
}

static inline void
find_idx_off_1st_rel(creal xi[3], int lg[3], creal og[3], creal shift,
		     creal dxi[3])
{
  for (int d = 0; d < 3; d++) {
    creal pos = xi[d] * dxi[d] + shift;
    lg[d] = fint(pos);
    og[d] = pos - lg[d];
  }
}

static inline void
find_idx_off_pos_1st(creal xi[3], int lg[3], creal og[3], creal pos[3], creal shift,
		     double xb[3], creal dxi[3])
{
  for (int d = 0; d < 3; d++) {
    pos[d] = (xi[d] - xb[d]) * dxi[d] + shift;
    lg[d] = fint(pos[d]);
    og[d] = pos[d] - lg[d];
  }
}

static inline void
find_idx_off_pos_1st_rel(creal xi[3], int lg[3], creal og[3], creal pos[3], creal shift,
			 creal dxi[3])
{
  for (int d = 0; d < 3; d++) {
    pos[d] = xi[d] * dxi[d] + shift;
    lg[d] = fint(pos[d]);
    og[d] = pos[d] - lg[d];
  }
}

#define INTERPOLATE_SETUP_1ST			\
  creal g0y = 1.f - og[1];			\
  creal g0z = 1.f - og[2];			\
  creal g1y = og[1];				\
  creal g1z = og[2];				\
						\
  creal h0y = 1.f - oh[1];			\
  creal h0z = 1.f - oh[2];			\
  creal h1y = oh[1];				\
  creal h1z = oh[2]


#define INTERPOLATE_FIELD_1ST(m, gy, gz)				\
    (gz##0z*(gy##0y*F3(pf, m, 0,l##gy[1]  ,l##gz[2]  ) +			\
	     gy##1y*F3(pf, m, 0,l##gy[1]+1,l##gz[2]  )) +			\
     gz##1z*(gy##0y*F3(pf, m, 0,l##gy[1]  ,l##gz[2]+1) +			\
	     gy##1y*F3(pf, m, 0,l##gy[1]+1,l##gz[2]+1)))

#ifdef F3_CURR

// ======================================================================
// current 1vb (yz)

static inline void
calc_dx1(creal dx1[2], creal x[2], creal dx[2], int off[2])
{
  if (off[1] == 0) {
    dx1[0] = .5f * off[0] - x[0];
    dx1[1] = dx[1] / dx[0] * dx1[0];
  } else {
    dx1[1] = .5f * off[1] - x[1];
    dx1[0] = dx[0] / dx[1] * dx1[1];
  }
}

static inline void
curr_2d_vb_cell(fields_curr_t *pf, int i[2], creal x[2], creal dx[2], creal fnq[2],
		creal dxt[2], int off[2])
{
  F3_CURR(pf, JYI, 0,i[0],i[1]  ) += fnq[0] * dx[0] * (.5f - x[1] - .5f * dx[1]);
  F3_CURR(pf, JYI, 0,i[0],i[1]+1) += fnq[0] * dx[0] * (.5f + x[1] + .5f * dx[1]);
  F3_CURR(pf, JZI, 0,i[0],i[1]  ) += fnq[1] * dx[1] * (.5f - x[0] - .5f * dx[0]);
  F3_CURR(pf, JZI, 0,i[0]+1,i[1]) += fnq[1] * dx[1] * (.5f + x[0] + .5f * dx[0]);
  if (dxt) {
    dxt[0] -= dx[0];
    dxt[1] -= dx[1];
    x[0] += dx[0] - off[0];
    x[1] += dx[1] - off[1];
    i[0] += off[0];
    i[1] += off[1];
  }
}

#endif

static void __unused
cache_fields_from_em(int p, fields_single_t *fld, fields_t *pf)
{
  struct psc_patch *patch = ppsc->patch + p;

  // FIXME, can do -1 .. 1?
  int ib[3] = { 0, -2, -2 };
  int ie[3] = { 1, patch->ldims[1] + 2, patch->ldims[2] + 2 };
  fields_single_alloc(fld, ib, ie, 9, 0); // JXI .. HZ
  for (int iz = -2; iz < patch->ldims[2] + 2; iz++) {
    for (int iy = -2; iy < patch->ldims[1] + 2; iy++) {
      F3_S(fld, EX, 0,iy,iz) = F3(pf, EX, 0,iy,iz);
      F3_S(fld, EY, 0,iy,iz) = F3(pf, EY, 0,iy,iz);
      F3_S(fld, EZ, 0,iy,iz) = F3(pf, EZ, 0,iy,iz);
      F3_S(fld, HX, 0,iy,iz) = F3(pf, HX, 0,iy,iz);
      F3_S(fld, HY, 0,iy,iz) = F3(pf, HY, 0,iy,iz);
      F3_S(fld, HZ, 0,iy,iz) = F3(pf, HZ, 0,iy,iz);
    }
  }
}

static void __unused
cache_fields_to_j(int p, fields_single_t *fld, fields_t *pf)
{
  struct psc_patch *patch = ppsc->patch + p;

  for (int iz = -2; iz < patch->ldims[2] + 2; iz++) {
    for (int iy = -2; iy < patch->ldims[1] + 2; iy++) {
      F3(pf, JXI, 0,iy,iz) += F3_S(fld, JXI, 0,iy,iz);
      F3(pf, JYI, 0,iy,iz) += F3_S(fld, JYI, 0,iy,iz);
      F3(pf, JZI, 0,iy,iz) += F3_S(fld, JZI, 0,iy,iz);
    }
  }
  fields_single_free(fld);
}

