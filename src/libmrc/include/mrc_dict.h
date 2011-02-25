
#ifndef MRC_DICT_H
#define MRC_DICT_H

#include <mrc_obj.h>

union param_u;

extern struct mrc_class mrc_class_mrc_dict;

MRC_OBJ_DEFINE_STANDARD_METHODS(mrc_dict, struct mrc_dict);
void mrc_dict_add_int(struct mrc_dict *dict, const char *name, int val);
void mrc_dict_add_bool(struct mrc_dict *dict, const char *name, bool val);
void mrc_dict_add_float(struct mrc_dict *dict, const char *name, float val);
void mrc_dict_add_double(struct mrc_dict *dict, const char *name, double val);
void mrc_dict_add_string(struct mrc_dict *dict, const char *name, const char *s);
void mrc_dict_add(struct mrc_dict *dict, int type, const char *name,
		  union param_u *pv);

#endif
