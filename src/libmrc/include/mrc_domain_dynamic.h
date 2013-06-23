//Forward declaration
struct mrc_domain_dynamic;
struct mrc_domain;

//
//void advance_grid_in_z(struct mrc_domain_dynamic* dynamic, struct mrc_domain* domain, int direction);

///@return gpatch of the newly created patch
int mrc_domain_dynamic_register_new_patch(struct mrc_domain* domain, int idx[3]);