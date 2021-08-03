PATH=/mnt/lustre/germaschewski/kaig/src/spack/opt/spack/linux-rhel7-zen/gcc-8.2.0/gcc-10.2.0-isfcn6dvyw2zpm75iiormpuwjsqxecnx/bin:/mnt/lustre/germaschewski/kaig/src/spack/opt/spack/cray-rhel7-zen2/gcc-10.2.0/git-2.29.0-uyrvs4qcqzy73qwema7xicg2gv23xwae/bin:/mnt/lustre/germaschewski/kaig/src/spack/opt/spack/cray-rhel7-zen2/gcc-10.2.0/cmake-3.18.4-uqgqfzos6p67ci37zz6zg645jtwunhhs/bin:/mnt/lustre/germaschewski/kaig/src/spack/opt/spack/cray-rhel7-zen2/gcc-10.2.0/openmpi-4.1.0-4majqkjuf6ppvo76eiewc3pfk44yo3tz/bin:$PATH
export CMAKE_PREFIX_PATH=/mnt/lustre/germaschewski/kaig/src/spack/opt/spack/cray-rhel7-zen2/gcc-10.2.0/openmpi-4.1.0-4majqkjuf6ppvo76eiewc3pfk44yo3tz:/mnt/lustre/germaschewski/kaig/src/spack/opt/spack/cray-rhel7-zen2/gcc-10.2.0/hdf5-1.10.7-jdmkhmdg3wgmolgqdpyfjwcmnfbut2sn
export CC=/mnt/lustre/germaschewski/kaig/src/spack/opt/spack/linux-rhel7-zen/gcc-8.2.0/gcc-10.2.0-isfcn6dvyw2zpm75iiormpuwjsqxecnx/bin/gcc
export CXX=/mnt/lustre/germaschewski/kaig/src/spack/opt/spack/linux-rhel7-zen/gcc-8.2.0/gcc-10.2.0-isfcn6dvyw2zpm75iiormpuwjsqxecnx/bin/g++

cmake -S .. -DCMAKE_BUILD_TYPE=Debug
