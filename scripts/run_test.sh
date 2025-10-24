#!/bin/bash

# Usage: scripts/run_test.sh name_of_test
# Example: scripts/run_test.sh test_injector_boundary_inflow
# Assumes you are in the root psc directory and have a build directory.

set -e
cd build
make $1
mkdir -p runs
cp ../bits/adios2cfg.xml runs/
cd runs
../src/libpsc/tests/$1 "${@:2}"
