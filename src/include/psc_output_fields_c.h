
#ifndef PSC_OUTPUT_FIELDS_C_H
#define PSC_OUTPUT_FIELDS_C_H

#include "psc_output_fields_private.h"
#include "psc_output_fields_item.h"

#include "fields3d.hxx"
#include "fields_item.hxx"

// ----------------------------------------------------------------------

enum {
  IO_TYPE_PFD,
  IO_TYPE_TFD,
  NR_IO_TYPES,
};

struct OutputFieldsItem
{
  OutputFieldsItem(PscFieldsItemBase item, const std::string& name,
		   std::vector<std::string>& comp_names, MfieldsBase& pfd,
		   MfieldsBase& tfd)
    : item(item), name(name), comp_names(comp_names), pfd(pfd), tfd(tfd)
  {}

  PscFieldsItemBase item;
  MfieldsBase& pfd;
  MfieldsBase& tfd;
  std::string name;
  std::vector<std::string> comp_names;
};

struct psc_output_fields_c {
  char *data_dir;
  char *output_fields;
  char *pfd_s;
  char *tfd_s;
  bool dowrite_pfield, dowrite_tfield;
  int pfield_first, tfield_first;
  int pfield_step, tfield_step;
  int tfield_length;
  int tfield_every;
  int rn[3];
  int rx[3];
	
  int pfield_next, tfield_next;
  // storage for output
  unsigned int naccum;
  std::vector<OutputFieldsItem> items;
  struct mrc_io *ios[NR_IO_TYPES];
};

#endif
