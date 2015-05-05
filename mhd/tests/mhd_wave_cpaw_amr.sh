#!/usr/bin/env bash
#
#PBS -l nodes=1:ppn=32
#PBS -l walltime=48:00:00
#PBS -j oe
#
# Arguments: ./mhd_wave_cpaw.sh [pdim [mode_num [ncells [v_par]]]]
#     mode_num (int): mode number in all dimensions of the box
#     ncells (int): number of grid cells in a patch (one dimension)
#     v_par (float): parallel flow, use -1.0 for a stationary wave
#
# Run Notes: Circularly polarized alfven wave. Use pdim and mode_num to
#            customize the domain, direction of wave propogation, and
#            wavelength
#

cd "${PBS_O_WORKDIR:-./}"

nprocs=2

mode=${1:-1}
ncells=${2:-16}
v_par=${3:-0.0}

MPIRUN="${MPIRUN:-mpirun}"
MPI_OPTS="-n ${nprocs} ${MPI_OPTS}"
GDB="${GDB:-gdb}"

run_name="cpaw_amr"

bin="./mhd_wave"
cmd="${MPIRUN} ${MPI_OPTS} ${bin}"
# cmd="${MPIRUN} ${MPI_OPTS} xterm -e ${GDB} --args ${bin}"

${cmd}                                                                       \
  --cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc \
  --ccccccccccccccccccccc  "GENERAL"                                         \
  --run                        ${run_name}                                   \
  --n_mhd_procs                ${nprocs}                                     \
  --cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc \
  --ccccccccccccccccccccc  "GRID"                                            \
  --mrc_domain_type            amr                                           \
  --mrc_crds_type              amr_uniform                                   \
  --crds_gen_x_type            uniform                                       \
  --crds_gen_y_type            uniform                                       \
  --crds_gen_z_type            uniform                                       \
  --mrc_domain_curve_type      bydim                                         \
  --ggcm_mhd_amr               4                                             \
  --mrc_domain_mx ${ncells} --mrc_domain_my ${ncells} --mrc_domain_mz  1     \
  --mrc_crds_lx    0.0      --mrc_crds_ly    0.0      --mrc_crds_lz    0.0   \
  --mrc_crds_hx    1.0      --mrc_crds_hy    1.0      --mrc_crds_hz    0.1   \
  --ggcm_mhd_crds_type         c                                             \
  --ggcm_mhd_crds_gen_type     mrc                                           \
  --cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc \
  --ccccccccccccccccccccc  "INITIAL CONDITION"                               \
  --ggcm_mhd_ic_type           cpaw                                          \
  --ggcm_mhd_ic_v_par          ${v_par}                                      \
  --ggcm_mhd_ic_polarization   1.0                                           \
  --ggcm_mhd_ic_mx ${mode} --ggcm_mhd_ic_my ${mode} --ggcm_mhd_ic_mz 0       \
  --cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc \
  --ccccccccccccccccccccc  "NUMERICS"                                        \
  --ggcm_mhd_step_type         c3_double                                     \
  --ggcm_mhd_do_vgrupd         0                                             \
  --cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc \
  --ccccccccccccccccccccc  "TIMESTEP"                                        \
  --ggcm_mhd_step_do_nwst      1                                             \
  --ggcm_mhd_thx               0.4                                           \
  --dtmin                      5e-5                                          \
  --mrc_ts_max_time            0.1                                           \
  --cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc \
  --ccccccccccccccccccccc  "SIM PARAMETERS"                                  \
  --ggcm_mhd_d_i               0.0                                           \
  --rrmin                      0.0                                           \
  --cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc \
  --ccccccccccccccccccccc  "OUTPUT"                                          \
  --mrc_io_type                xdmf2                                         \
  --ggcm_mhd_step_debug_dump   0                                             \
  --mrc_ts_output_every_time   0.05                                          \
  --mrc_io_sw                  0                                             \
  --ggcm_mhd_diag_fields       rr:pp:v:b:b1:j:e_cc:divb                      \
  --cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc \
  --ccccccccccccccccccccc  "MISC"                                            \
  --ggcm_mhd_do_badval_checks  1                                             \
  --monitor_conservation       0                                             \
  --cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc \
  2>&1 | tee ${run_name}.log

errcode=$?

# plot quick and dirty time series with:
# viscid_ts --slice x=0.0,y=0.0,z=0.0 -p vx -p vy -p vz -p bx -p by -p bz \
#           cpaw_amr.3d.xdmf

# plot final time step:
# viscid_2d --slice z=0.0 -o p_y -p vx -p vy -p vz -p bx -p by -p bz \
#           -t -1 cpaw_amr.3d.xdmf

# plot movie:
# viscid_2d --slice z=0.0 -o p_y -p vx -p vy -p vz -p bx -p by -p bz \
#           -a cpaw_amr.mp4 --np 2 cpaw_amr.3d.xdmf

exit ${errcode}

##
## EOF
##
