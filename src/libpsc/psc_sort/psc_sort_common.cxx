
using real_t = mparticles_t::real_t;

// ======================================================================
// quicksort

static mparticles_t::patch_t *s_prts;

static int
compare(const void *_a, const void *_b)
{
  const particle_t *a = (const particle_t *) _a, *b = (const particle_t *) _b;

  int cni_a = s_prts->validCellIndex(*a);
  int cni_b = s_prts->validCellIndex(*b);
  if (cni_a < cni_b) {
    return -1;
  } else if (cni_a == cni_b) {
    return 0;
  } else {
    return 1;
  }
}

static void
psc_sort_qsort_run(struct psc_sort *sort, struct psc_mparticles *mprts_base)
{
  // FIXME, using a C++ sort would be much nicer...
  static int pr;
  if (!pr) {
    pr = prof_register("qsort_sort", 1., 0, 0);
  }
  mparticles_t mprts = mprts_base->get_as<mparticles_t>();

  prof_start(pr);
  for (int p = 0; p < mprts->n_patches(); p++) {
    s_prts = &mprts[p];
    qsort(&*s_prts->begin(), s_prts->size(), sizeof(*s_prts->begin()), compare);
  }
  prof_stop(pr);

  mprts.put_as(mprts_base);
}

// ======================================================================
// counting sort

static void
psc_sort_countsort_run(struct psc_sort *sort, struct psc_mparticles *mprts_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("countsort_sort", 1., 0, 0);
  }

  mparticles_t mprts = mprts_base->get_as<mparticles_t>();

  prof_start(pr);
  for (int p = 0; p < mprts->n_patches(); p++) {
    mparticles_t::patch_t& prts = mprts[p];
    unsigned int n_prts = prts.size();

    unsigned int n_cells = prts.pi_.n_cells_;
    unsigned int *cnts = (unsigned int *) malloc(n_cells * sizeof(*cnts));
    memset(cnts, 0, n_cells * sizeof(*cnts));
    
    // count
    PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
      int cni = prts.validCellIndex(*prt_iter);
      cnts[cni]++;
    }
    
    // calc offsets
    int cur = 0;
    for (int i = 0; i < n_cells; i++) {
      int n = cnts[i];
      cnts[i] = cur;
    cur += n;
    }
    assert(cur == n_prts);
    
    // move into new position
    particle_t *particles2 = (particle_t *) malloc(n_prts * sizeof(*particles2));
    PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
      unsigned int cni = prts.validCellIndex(*prt_iter);
      memcpy(&particles2[cnts[cni]], &*prt_iter, sizeof(*particles2));
      cnts[cni]++;
    }
    
    // back to in-place
    memcpy(&*prts.begin(), particles2, n_prts * sizeof(*particles2));
    
    free(particles2);
    free(cnts);
  }
  
  prof_stop(pr);

  mprts.put_as(mprts_base);
}


// ======================================================================
// counting sort 2
// use a separate array of cell indices

template<typename M>
struct psc_sort_countsort2
{
  using mparticles_t = M;
  
  void operator()(mparticles_t mprts)
  {
    for (int p = 0; p < mprts->n_patches(); p++) {
      auto& prts = mprts[p];
      unsigned int n_prts = prts.size();
      
      unsigned int n_cells = prts.pi_.n_cells_;
      unsigned int *cnis = (unsigned int *) malloc(n_prts * sizeof(*cnis));
      // FIXME, might as well merge counting here, too
      int i = 0;
      for (auto prt_iter = prts.begin(); prt_iter != prts.end(); ++prt_iter, ++i) {
	int cni = prts.validCellIndex(*prt_iter);
	assert(cni >= 0);
	cnis[i] = cni;
      }
      
      unsigned int *cnts = (unsigned int *) malloc(n_cells * sizeof(*cnts));
      memset(cnts, 0, n_cells * sizeof(*cnts));
	
      // count
      for (int i = 0; i < n_prts; i++) {
	unsigned int cni = cnis[i];
	cnts[cni]++;
      }
	
      // calc offsets
      int cur = 0;
      for (int i = 0; i < n_cells; i++) {
	int n = cnts[i];
	cnts[i] = cur;
	cur += n;
      }
      assert(cur == n_prts);
      
      // move into new position
      particle_t *particles2 = (particle_t *) malloc(n_prts * sizeof(*particles2));
      for (int i = 0; i < n_prts; i++) {
	unsigned int cni = cnis[i];
	int n = 1;
	while (i+n < n_prts && cnis[i+n] == cni) {
	  n++;
	}
	memcpy(&particles2[cnts[cni]], &prts[i], n * sizeof(*particles2));
	cnts[cni] += n;
	i += n - 1;
      }
      
      // back to in-place
      memcpy(&*prts.begin(), particles2, n_prts * sizeof(*particles2));
      
      free(particles2);
      free(cnis);
      free(cnts);
    }
  }
 
  static void run(struct psc_sort *sort, struct psc_mparticles *mprts_base)
  {
    static int pr;
    if (!pr) {
      pr = prof_register("countsort2_sort", 1., 0, 0);
    }
    
    mparticles_t mprts = mprts_base->get_as<mparticles_t>();
    
    prof_start(pr);
    psc_sort_countsort2& sub = *mrc_to_subobj(sort, psc_sort_countsort2);
    sub(mprts);
    prof_stop(pr);
    
    mprts.put_as(mprts_base);
  }
};

