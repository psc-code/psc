# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class Rmm(CMakePackage):
    """RMM: RAPIDS Memory Manager"""

    homepage = "https://github.com/rapidsai/rmm"
    git      = "https://github.com/rapidsai/rmm"

    maintainers = ['germasch']

    version('0.16.0', tag='v0.16.0', preferred=True)

    variant('tests', default=False, description="Build tests")
    
    depends_on('cuda')
    depends_on('spdlog@1.7.0')
    #depends_on('thrust@1.10.0:')
    depends_on('googletest@1.10.0 +gmock', type='build', when='+tests')

    patch('rmm-summit.patch', when='@0.16.0:')

    def cmake_args(self):
        spec = self.spec

        args = []
        args += ['-DBUILD_TESTS:BOOL={}'.format(
            'ON' if '+tests' in spec else 'OFF')]

        return args
    
