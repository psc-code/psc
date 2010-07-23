void
init_output_fields(struct psc_extra_fields *f);

void
free_output_fields(struct psc_extra_fields *f);

void
calculate_pfields(struct psc_extra_fields *p);

void
accumulate_tfields(struct psc_extra_fields *p, struct psc_extra_fields *t);

void
reset_fields(struct psc_extra_fields *f);

// convert accumulated values to correct temporal mean
// (divide by naccum)
void
mean_tfields(struct psc_extra_fields *f);
