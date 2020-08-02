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

    depends_on('hdf5@1.8.0: +hl')

