# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class Thrust(CMakePackage):
    """Thrust is a parallel algorithms library
    which resembles the C++ Standard Template Library (STL)."""

    homepage = "https://thrust.github.io"
    url      = "https://github.com/NVIDIA/thrust/archive/1.10.0.tar.gz"
    git      = "https://github.com/NVIDIA/thrust"

    version('develop', branch='main', submodules=True)
    # FIXME, releases should be based on tarball, but need cub then
    version('1.12.0', tag='1.12.0', submodules=True)
    version('1.10.0', tag='1.10.0', submodules=True)
    #version('1.12.0', sha256='d68f89937ed1a0dadafd4f77d3e65ab9a5baa48dd5fa5698f8d3f7e3d6749da2')
    #version('1.11.0', sha256='c65211a66fe0dce9f5470cc8fcb80ae6b04da12fd94707ed63ddddbd82fa9444')
    version('1.9.10', sha256='5071c995a03e97e2bcfe0c5cc265330160316733336bb87dfcea3fc381065316')
    version('1.9.9', sha256='74740b94e0c62e1e54f9880cf1eeb2d0bfb2203ae35fd72ece608f9b8e073ef7')
    version('1.9.8', sha256='d014396a2128417bd1421ba80d2601422577094c0ff727c60bd3c6eb4856af9b')
    version('1.9.7', sha256='72519f7f6b5d28953e5086253bbcf5b10decde913ddeb93872dc51522bdfad2b')
    version('1.9.6', sha256='67e937c31d279cec9ad2c54a4f66e479dfbe506ceb0253f611a54323fb4c1cfb')
    version('1.9.5', sha256='d155dc2a260fe0c75c63c185fa4c4b4c6c5b7c444fcdac7109bb71941c9603f1')
    version('1.9.4', sha256='41931a7d73331fc39c6bea56d1eb8d4d8bbf7c73688979bbdab0e55772f538d1')
    version('1.9.3', sha256='92482ad0219cd2d727766f42a4fc952d7c5fd0183c5e201d9a117568387b4fd1')
    version('1.9.2', sha256='1fb1272be9e8c28973f5c39eb230d1914375ef38bcaacf09a3fa51c6b710b756')
    version('1.9.1', sha256='7cf59bf42a7b05bc6799c88269bf41eb637ca2897726a5ade334a1b8b4579ef1')
    version('1.9.0', sha256='a98cf59fc145dd161471291d4816f399b809eb0db2f7085acc7e3ebc06558b37')
    version('1.8.2', sha256='83bc9e7b769daa04324c986eeaf48fcb53c2dda26bcc77cb3c07f4b1c359feb8')

    patch('1323-mod.patch', when='@develop')
    patch('thrust-summit.patch')
    
    patch('https://github.com/NVIDIA/thrust/pull/1297.patch', when='@1.10.0',
              sha256='30f8c67a0e8fd8584a233dc95e2ff68a37cfc58a55872381bae7281c8c74045a')
    patch('https://github.com/NVIDIA/thrust/pull/1298.patch', when='@1.10.0',
              sha256='0f27f6883e0f16d2d83f8816ab09dad4ad890aed8556ba1c1227e5339ed115c0')
    patch('https://github.com/NVIDIA/thrust/pull/1300.patch', when='@1.10.0',
              sha256='dbcbe70de701c5f64c9307cc921c2de06924843ef187fe59f69d11e41dda5211')

    depends_on('cuda')
    
    def build(self, spec, prefix):
        pass
