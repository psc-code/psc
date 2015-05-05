#!/usr/bin/env bash

has_viscid=$(python -c 'import viscid' &> /dev/null; echo ${?})

./mhd_wave_cpaw.sh x 1 64 -1.0 &> /dev/null
errcode=$?

if [ ${has_viscid} -eq 0 ]; then
  python check_wave_cpaw.py cpaw_x_64.3d.xdmf
  errcode=$?
else
  echo "Viscid python module not found, not checking result." >&2
fi

if [ $errcode -eq 0 ]; then
  rm cpaw_x_64.3d.*xdmf
  rm cpaw_x_64.3d.*h5
fi

exit ${errcode}

##
## EOF
##
