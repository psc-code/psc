
#ifndef VPIC_IFACE_H
#define VPIC_IFACE_H

#ifdef __cplusplus
extern "C" {
#endif
#if 0 // hack to fix indentation
}
#endif

struct vpic_info {
  double dx, dy, dz;
  double dt;
  double c;
  double eps0;
};

void vpic_base_init(struct vpic_info *info);
void vpic_base_integrate();

bool vpic_done();

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
