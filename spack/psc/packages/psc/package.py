# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class Psc(CMakePackage):
    """The PSC particle-in-cell code."""

    homepage = "https://psc.readthedocs.io/"
    url      = "https://psc.readthedocs.io/fake.tar.gz"
    git      = "https://github.com/psc-code/psc"

    maintainers = ['germasch']

    version('develop', branch='master')

    variant('adios2', default=True,
            description='Enable ADIOS2 support')
    variant('cuda', default=False,
            description='Enable CUDA support')
    variant('nvtx', default=False,
            description='Enable NVTX profiling support')
    variant('rmm', default=False,
            description='Enable RMM memory manager support')
    variant('tests', default=False,
            description='Build with unit testing')
    
    depends_on('cmake@3.17.0:')

    depends_on('hdf5 +hl')
    depends_on('adios2@2.4.0:', when='+adios2')
    depends_on('cuda', when='+cuda')
    depends_on('thrust@1.10.0:', when='+cuda')
    depends_on('cuda', when='+nvtx')
    depends_on('rmm', when='+rmm')

    def cmake_args(self):
        spec = self.spec
        
        args = []
        args += ['-DUSE_CUDA={}'.format(
            'ON' if '+cuda' in self.spec else 'OFF')]
        args += ['-DPSC_USE_NVTX={}'.format(
            'ON' if '+nvtx' in self.spec else 'OFF')]
        args += ['-DPSC_USE_RMM={}'.format(
            'ON' if '+rmm' in self.spec else 'OFF')]
        args += ['-DBUILD_TESTING={}'.format(
            'ON' if '+tests' in self.spec else 'OFF')]
        if 'thrust' in spec:
            args += ['-DThrust_DIR={}/include/thrust/cmake'.format(spec['thrust'].prefix)]

        return args
