#! /bin/sh

openmpirun -n 4 ./test_mrc_domain_multi --npx 3 --npy 2 --mrc_io_type xdmf2
