
#include "cuda_mfields.h"
#include "cuda_bits.h"
#include "cuda_base.cuh"

#include "fields.hxx"

#include <cstdio>
#include <cassert>

// ======================================================================
// cuda_mfields

// ----------------------------------------------------------------------
// to_json

mrc_json_t cuda_mfields::to_json()
{
  mrc_json_t json = mrc_json_object_new(9);
  mrc_json_object_push_integer(json, "n_patches", n_patches());
  mrc_json_object_push_integer(json, "n_fields", n_comps());

  // mrc_json_object_push(json, "ib", mrc_json_integer_array_new(3, ib()));
  // mrc_json_object_push(json, "im", mrc_json_integer_array_new(3, im()));

  mrc_json_t json_flds = mrc_json_object_new(2);
  mrc_json_object_push(json, "flds", json_flds);
  mrc_json_object_push_boolean(json_flds, "__field5d__", true);
  mrc_json_t json_flds_patches = mrc_json_array_new(n_patches());
  mrc_json_object_push(json_flds, "data", json_flds_patches);

  auto h_mflds = hostMirror(*this);
  copy(*this, h_mflds);
  for (int p = 0; p < n_patches(); p++) {
    auto flds = make_Fields3d<dim_xyz>(h_mflds[p]);

    mrc_json_t json_flds_comps = mrc_json_array_new(n_comps());
    mrc_json_array_push(json_flds_patches, json_flds_comps);
    for (int m = 0; m < n_comps(); m++) {
      mrc_json_t json_fld_z = mrc_json_array_new(im(2));
      mrc_json_array_push(json_flds_comps, json_fld_z);
      for (int k = ib(2); k < ib(2) + im(2); k++) {
        mrc_json_t json_fld_y = mrc_json_array_new(im(1));
        mrc_json_array_push(json_fld_z, json_fld_y);
        for (int j = ib(1); j < ib(1) + im(1); j++) {
          mrc_json_t json_fld_x = mrc_json_array_new(im(0));
          mrc_json_array_push(json_fld_y, json_fld_x);
          for (int i = ib(0); i < ib(0) + im(0); i++) {
            mrc_json_array_push_double(json_fld_x, flds(m, i, j, k));
          }
        }
      }
    }
  }

  return json;
}

// ----------------------------------------------------------------------
// dump

void cuda_mfields::dump(const char* filename)
{
  mrc_json_t json = to_json();

  const char* buf = mrc_json_to_string(json);
  if (filename) {
    FILE* file = fopen(filename, "w");
    assert(file);
    fwrite(buf, 1, strlen(buf), file);
    fclose(file);
  } else {
    printf("cuda_mfields (json):\n%s\n", buf);
  }
  free((void*)buf);

  // FIXME free json
}

// ----------------------------------------------------------------------
// cast to DMFields

cuda_mfields::operator DMFields()
{
  return DMFields{box(), n_comps(), n_patches(), storage().data().get()};
}

// ----------------------------------------------------------------------
// operator[]

DFields cuda_mfields::operator[](int p) const
{
  return static_cast<DMFields>(const_cast<cuda_mfields&>(*this))[p];
}

HMFields hostMirror(cuda_mfields& cmflds)
{
  return HMFields{cmflds.grid(), cmflds.n_comps(), -cmflds.ib()};
}

HMFields hostMirror(const cuda_mfields& cmflds)
{
  return HMFields{cmflds.grid(), cmflds.n_comps(), -cmflds.ib()};
}

void copy(const cuda_mfields& cmflds, HMFields& hmflds)
{
  static int pr;
  if (!pr) {
    pr = prof_register("cmflds to host", 1., 0, 0);
  }
  prof_start(pr);
  thrust::copy(cmflds.storage().data(),
               cmflds.storage().data() + cmflds.storage().size(),
               hmflds.storage().data());
  prof_stop(pr);
}

void copy(const HMFields& hmflds, cuda_mfields& cmflds)
{
  static int pr;
  if (!pr) {
    pr = prof_register("cmflds from host", 1., 0, 0);
  }
  prof_start(pr);
  thrust::copy(hmflds.storage().data(),
               hmflds.storage().data() + hmflds.storage().size(),
               cmflds.storage().data());
  prof_stop(pr);
}
