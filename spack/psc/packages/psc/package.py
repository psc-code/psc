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

    depends_on('hdf5@1.8.0:1.8.999 +hl')
    depends_on('adios2@2.4.0:', when='+adios2')

