
// OPT, this could probably be optimized by using a struct, and doing one
// copy

__constant__ static real d_dt, d_dxi[3], d_dqs, d_fnqs, d_fnqys, d_fnqzs;
__constant__ static int d_mx[3], d_ilg[3], d_ilo[3];
__constant__ static int d_b_mx[3];

EXTERN_C void
PFX(set_constants)(particles_cuda_t *pp, fields_cuda_t *pf)
{
  real __dt = ppsc->dt, __dqs = .5f * ppsc->coeff.eta * ppsc->dt;
  real fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  real fnqys = ppsc->dx[1] * fnqs / ppsc->dt;
  real fnqzs = ppsc->dx[2] * fnqs / ppsc->dt;
  check(cudaMemcpyToSymbol(d_dt, &__dt, sizeof(d_dt)));
  check(cudaMemcpyToSymbol(d_dqs, &__dqs, sizeof(d_dqs)));
  check(cudaMemcpyToSymbol(d_fnqs, &fnqs, sizeof(d_fnqs)));
  check(cudaMemcpyToSymbol(d_fnqys, &fnqys, sizeof(d_fnqys)));
  check(cudaMemcpyToSymbol(d_fnqzs, &fnqzs, sizeof(d_fnqzs)));
  real __dxi[3] = { 1.f / ppsc->dx[0], 1.f / ppsc->dx[1], 1.f / ppsc->dx[2] };
  check(cudaMemcpyToSymbol(d_dxi, __dxi, sizeof(d_dxi)));
  check(cudaMemcpyToSymbol(d_mx, pf->im, sizeof(d_mx)));
  check(cudaMemcpyToSymbol(d_ilg, pf->ib, sizeof(d_ilg)));
  real ilo[3] = { pf->ib[0] + ppsc->ibn[0], pf->ib[1] + ppsc->ibn[1], pf->ib[2] + ppsc->ibn[2] };
  check(cudaMemcpyToSymbol(d_ilo, ilo, sizeof(d_ilo)));
  check(cudaMemcpyToSymbol(d_b_mx, pp->b_mx, sizeof(d_mx)));
}

