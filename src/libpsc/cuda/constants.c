
struct d_consts {
  real dt;
  real dxi[3];
  real dqs;
  real fnqs;
  real fnqys, fnqzs;
  int mx[3];
  int ilg[3];
  int ilo[3];
  //  int b_mx[3];
};

__constant__ static struct d_consts d_consts;

EXTERN_C void
PFX(set_constants)(struct psc_mfields *mflds, int p)
{
  struct cuda_mfields *cmflds = psc_mfields_cuda(mflds)->cmflds;
  
  struct d_consts consts;
  consts.dt     = ppsc->dt;
  consts.dqs    = .5f * ppsc->coeff.eta * ppsc->dt;
  consts.fnqs   = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  consts.fnqys  = ppsc->patch[p].dx[1] * consts.fnqs / ppsc->dt;
  consts.fnqzs  = ppsc->patch[p].dx[2] * consts.fnqs / ppsc->dt;
  for (int d = 0; d < 3; d++) {
    consts.dxi[d] = 1.f / ppsc->patch[p].dx[d];
    consts.mx[d]  = cmflds->im[d];
    consts.ilg[d] = cmflds->ib[d];
    consts.ilo[d] = cmflds->ib[d] + ppsc->ibn[d];
  }

#if 0
  if (prts) {
    struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(prts->mprts);
    for (int d = 0; d < 3; d++) {
      consts.b_mx[d] = mprts_cuda->b_mx[d];
    }
  }
#endif

  cudaError_t ierr;
  ierr = cudaMemcpyToSymbol(d_consts, &consts, sizeof(consts)); cudaCheck(ierr);
}

