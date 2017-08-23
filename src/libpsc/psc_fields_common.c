
// ----------------------------------------------------------------------
// psc_fields_setup

static void
PFX(setup)(struct psc_fields *pf)
{
  unsigned int size = 1;
  for (int d = 0; d < 3; d++) {
    size *= pf->im[d];
  }

#if PSC_FIELDS_AS_FORTRAN
  fields_real_t **flds = calloc(pf->nr_comp, sizeof(*flds));
  flds[0] = calloc(size * pf->nr_comp, sizeof(flds[0]));
  for (int i = 1; i < pf->nr_comp; i++) {
    flds[i] = flds[0] + i * size;
  }
  pf->data = flds;
#elif PSC_FIELDS_AS_C && defined(USE_CBE)
  // The Cell processor translation can use the C fields with one modification:
  // the data needs to be 128 byte aligned (to speed off-loading to spes). This
  // change is roughly put in below.
  void *m;
  int ierr = posix_memalign(&m, 128, nr_comp * size * sizeof(fields_real_t));
  assert(ierr == 0);
  pf->flds = m; 
#else
  pf->data = calloc(pf->nr_comp * size, sizeof(fields_real_t));
#endif
}

// ----------------------------------------------------------------------
// psc_fields_destroy

static void
PFX(destroy)(struct psc_fields *pf)
{
#if PSC_FIELDS_AS_FORTRAN
  fields_real_t **flds = pf->data;
  free(flds[0]);

  for (int i = 0; i < pf->nr_comp; i++) {
    flds[i] = NULL;
  }
  free(flds);
#else
  free(pf->data);
#endif
}

// ----------------------------------------------------------------------
// psc_fields_zero_comp

static void
PFX(zero_comp)(struct psc_fields *pf, int m)
{
  memset(&F3(pf, m, pf->ib[0], pf->ib[1], pf->ib[2]), 0,
	 pf->im[0] * pf->im[1] * pf->im[2] * sizeof(fields_real_t));
}

// ----------------------------------------------------------------------
// psc_fields_set_comp

static void
PFX(set_comp)(struct psc_fields *pf, int m, double _val)
{
  fields_real_t val = _val;

  for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
    for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
      for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	F3(pf, m, jx, jy, jz) = val;
      }
    }
  }
}

