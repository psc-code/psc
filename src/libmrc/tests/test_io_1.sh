#! /bin/sh
openmpirun -n 2 ./test_io --npx 2 --mrc_io_type xdmf_to_one
