
#ifndef PSC_BALANCE_H
#define PSC_BALANCE_H

#include <mrc_obj.h>

#include "psc.h"

BEGIN_C_DECLS

MRC_CLASS_DECLARE(psc_balance, struct psc_balance);

void psc_balance_initial(struct psc_balance *bal, struct psc *psc,
			 uint*& n_prts_by_patch);
void psc_balance_run(struct psc_balance *bal, struct psc *psc);

extern int psc_balance_generation_cnt;

END_C_DECLS;

#endif
