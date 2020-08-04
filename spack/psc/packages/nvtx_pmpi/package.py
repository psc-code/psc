# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class NvtxPmpi(CMakePackage):
    """FIXME: Put a proper description of your package here."""

    # FIXME: Add a proper url for your package's homepage here.
    homepage = "https://www.example.com"
    url      = "nvtx_pmpi"
    git      = "https://github.com/germasch/cuda-profiler.git"

    # FIXME: Add a list of GitHub accounts to
    # notify when the package is updated.
    maintainers = ['germasch']

    # FIXME: Add proper versions and checksums here.
    version('0.0.1', branch='pr/cmake')

    depends_on('mpi')
    depends_on('cuda')
    depends_on('python@:2.7.999', type='build')
    depends_on('cmake@3.17.0:'  , type='build')

    root_cmakelists_dir = 'nvtx_pmpi_wrappers'
