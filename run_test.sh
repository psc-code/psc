#!/bin/bash
set -e
cd build
make $1
mkdir -p runs
cp ../bits/adios2cfg.xml runs/
cd runs
../src/libpsc/tests/$1 "${@:2}"
