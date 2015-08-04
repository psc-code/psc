
#define MAX_NR_KINDS (10)

struct params_1vb {
  particle_real_t dt;
  particle_real_t fnqs, fnqxs, fnqys, fnqzs;
  particle_real_t dxi[3];
  particle_real_t dq_kind[MAX_NR_KINDS];
  particle_real_t fnqx_kind[MAX_NR_KINDS];
  particle_real_t fnqy_kind[MAX_NR_KINDS];
  particle_real_t fnqz_kind[MAX_NR_KINDS];
};

static struct params_1vb prm;

static void
params_1vb_set(struct psc *psc, int p)
{
  prm.dt = ppsc->dt;
  prm.fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
#ifdef VB_2D
#if !(DIM == DIM_YZ)
#error inc_params.c: VB_2D only works for DIM_YZ
#endif
  prm.fnqxs = prm.fnqs;
#else
  prm.fnqxs = ppsc->patch[p].dx[0] * prm.fnqs / prm.dt;
#endif
  prm.fnqys = ppsc->patch[p].dx[1] * prm.fnqs / prm.dt;
  prm.fnqzs = ppsc->patch[p].dx[2] * prm.fnqs / prm.dt;
  for (int d = 0; d < 3; d++) {
    prm.dxi[d] = 1.f / ppsc->patch[p].dx[d];
  }

  assert(ppsc->nr_kinds <= MAX_NR_KINDS);
  for (int k = 0; k < ppsc->nr_kinds; k++) {
    prm.dq_kind[k] = .5f * ppsc->coeff.eta * prm.dt * ppsc->kinds[k].q / ppsc->kinds[k].m;
    prm.fnqx_kind[k] = prm.fnqxs * ppsc->kinds[k].q;
    prm.fnqy_kind[k] = prm.fnqys * ppsc->kinds[k].q;
    prm.fnqz_kind[k] = prm.fnqzs * ppsc->kinds[k].q;
  }
}
