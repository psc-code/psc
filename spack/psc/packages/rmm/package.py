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

    version('cmake', branch='pr/cmake', git='https://github.com/germasch/rmm')

    depends_on('cuda')

