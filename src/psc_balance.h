
#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_balance, struct psc_balance);

void psc_balance_initial(struct psc_balance *bal, struct psc *psc,
			 int **p_nr_particles_by_patch);
void psc_balance_run(struct psc_balance *bal, struct psc *psc);
