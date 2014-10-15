#!/usr/bin/env bash
# This script uses the libmrc grid generation utility to generate a
# grid and if successful, interactively plots it.
#
# required:
#     - you must define LIBMRC_DIR to point to somewhere with a
#       compiled libmrc
#     - plots require python and a verison of matplotlib that has GridSpec.

if [ -z ${LIBMRC_DIR} ]; then
    if [ -z ${OPENGGCMDIR} ]; then
        LIBMRC_DIR=".."
    else
        LIBMRC_DIR="${OPENGGCMDIR}/target-build/csrc/libmrc"
    fi
fi

plot_bin="${LIBMRC_DIR}/util/plot_grid.py"
plot_opts=""
# plot_opts="${plot_opts} --sample ~/dev/stage/txsmall/target/*.3d.000001.xdmf"
# plot_opts="${plot_opts} -p J=sqrt(jx**2+jy**2+jz**2) -o cmap_Greens"

# ################################################
# #  generate a 128x64x64 grid with a tanh shape
# ${LIBMRC_DIR}/util/gen_mrc_crds          \
#   --mrc_crds_type     rectilinear        \
#                                          \
#   --crds_gen_x_type    ggcm_x_tanh       \
#   --mrc_domain_mx      128               \
#   --mrc_crds_lx       -30.0              \
#   --mrc_crds_hx        300.0             \
#   --crds_gen_x_x1     -26.0              \
#   --crds_gen_x_x3      10.0              \
#   --crds_gen_x_dmm     8.0               \
#   --crds_gen_x_x5      80.0              \
#   --crds_gen_x_h0      2.0               \
#   --crds_gen_x_hn      43.0              \
#   --crds_gen_x_hm      1.7               \
#   --crds_gen_x_hmm     0.7               \
#   --crds_gen_x_b1      0.15              \
#   --crds_gen_x_b2      0.025             \
#   --crds_gen_x_b3      0.3               \
#                                          \
#   --crds_gen_y_type            ggcm_yz   \
#   --mrc_domain_my              64        \
#   --mrc_crds_ly               -50.0      \
#   --mrc_crds_hy                50.0      \
#   --crds_gen_y_center_spacing  0.4       \
#   --crds_gen_y_center_shift    0.0       \
#   --crds_gen_y_xn              2.0       \
#   --crds_gen_y_xm              0.5       \
#                                          \
#   --crds_gen_z_type            ggcm_yz   \
#   --mrc_domain_mz              64        \
#   --mrc_crds_lz               -50.0      \
#   --mrc_crds_hz                50.0      \
#   --crds_gen_z_center_spacing  0.4       \
#   --crds_gen_z_center_shift    0.0       \
#   --crds_gen_z_xn              2.0       \
#   --crds_gen_z_xm              0.5
# if [ $? -eq 0 ]; then
#   python ${plot_bin} ${plot_opts} test.grid2 test.hgrid2
# fi

################################################
#  generate a 128x64x64 grid with a cubic shape
${LIBMRC_DIR}/util/gen_mrc_crds          \
  --mrc_crds_type     rectilinear        \
                                         \
  --crds_gen_x_type    ggcm_x_cubic      \
  --mrc_domain_mx      800               \
  --mrc_crds_lx       -30.0              \
  --mrc_crds_hx        300.0             \
  --crds_gen_x_w0      1.0               \
  --crds_gen_x_w1      150.0             \
  --crds_gen_x_a1      4.0               \
  --crds_gen_x_b1      400.0             \
  --crds_gen_x_w2      2.0               \
  --crds_gen_x_a2     -8.0               \
  --crds_gen_x_b2     -30.0              \
                                         \
  --crds_gen_y_type            ggcm_yz   \
  --mrc_domain_my              256       \
  --mrc_crds_ly               -50.0      \
  --mrc_crds_hy                50.0      \
  --crds_gen_y_center_spacing  0.05      \
  --crds_gen_y_center_shift    0.0       \
  --crds_gen_y_xn              2.0       \
  --crds_gen_y_xm              0.5       \
                                         \
  --crds_gen_z_type            ggcm_yz   \
  --mrc_domain_mz              256       \
  --mrc_crds_lz               -50.0      \
  --mrc_crds_hz                50.0      \
  --crds_gen_z_center_spacing  0.05      \
  --crds_gen_z_center_shift    0.0       \
  --crds_gen_z_xn              2.0       \
  --crds_gen_z_xm              0.5
if [ $? -eq 0 ]; then
  python ${plot_bin} ${plot_opts} test.grid2 test.hgrid2
fi

##
## EOF
##
