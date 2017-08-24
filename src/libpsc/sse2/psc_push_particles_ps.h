
#ifndef PSC_PUSH_PARTICLES_PS_H
#define PSC_PUSH_PARTICLES_PS_H

#include "psc.h"
#include "psc_particles_single.h"

#include "psc_fields_as_single.h"

void psc_push_particles_1vb_ps_push_a_yz(struct psc_push_particles *push,
					 struct psc_mparticles *mprts_base,
					 struct psc_mfields *mflds_base);

void psc_push_particles_1vb_ps2_push_a_yz(struct psc_push_particles *push,
					  struct psc_mparticles *mprts_base,
					  struct psc_mfields *mflds_base);

// ======================================================================

// FIXME, use regular psc_fields?

typedef struct {
  float f00, f01, f10, f11;
} fields_ip_real_t;

typedef struct {
  fields_ip_real_t *flds;
  int ib[3], im[3]; //> lower bounds and length per direction
  int nr_comp; //> nr of components
  int first_comp; // first component
} fields_ip_t;

#define F3_OFF_IP(pf, fldnr, jx,jy,jz)					\
  ((((((fldnr - (pf)->first_comp)					\
       * (pf)->im[2] + ((jz)-(pf)->ib[2]))				\
      * (pf)->im[1] + ((jy)-(pf)->ib[1]))				\
     * (pf)->im[0] + ((jx)-(pf)->ib[0]))))

#ifndef BOUNDS_CHECK

#define F3_IP(pf, fldnr, jx,jy,jz)		\
  ((pf)->flds[F3_OFF_IP(pf, fldnr, jx,jy,jz)])

#else

#define F3_IP(pf, fldnr, jx,jy,jz)					\
  (*({int off = F3_OFF_IP(pf, fldnr, jx,jy,jz);				\
      assert(fldnr >= (pf)->first_comp && fldnr < (pf)->first_comp + (pf)->nr_comp); \
      assert(jx >= (pf)->ib[0] && jx < (pf)->ib[0] + (pf)->im[0]);	\
      assert(jy >= (pf)->ib[1] && jy < (pf)->ib[1] + (pf)->im[1]);	\
      assert(jz >= (pf)->ib[2] && jz < (pf)->ib[2] + (pf)->im[2]);	\
      &((pf)->flds[off]);						\
    }))

#endif

void sb0_ps_1vb_yz_p(int p, fields_ip_t *pf, struct psc_particles *pps,
		     int n_start);
void sb2_ps_1vb_yz_p(int p, fields_ip_t *pf, struct psc_particles *pps,
		     int n_start);
void sb0_ps_1vb_yz_xx_jxyz(int p, fields_ip_t *pf, struct psc_particles *pps,
			   int n_start);
void sb2_ps_1vb_yz_xx_jxyz(int p, fields_ip_t *pf, struct psc_particles *pps,
			   int n_start);

void sb0_ps_1vb_yz_pxx_jxyz(int p, fields_ip_t *pf, struct psc_particles *pps,
			    int n_start);
void sb2_ps_1vb_yz_pxx_jxyz(int p, fields_ip_t *pf, struct psc_particles *pps,
			    int n_start);
void sb0_ps2_1vb_yz_pxx_jxyz(int p, fields_ip_t *pf, struct psc_particles *pps,
			     int n_start);
void sb2_ps2_1vb_yz_pxx_jxyz(int p, fields_ip_t *pf, struct psc_particles *pps,
			     int n_start);


#endif
