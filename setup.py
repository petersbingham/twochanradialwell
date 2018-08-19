# -*- coding: utf-8 -*-

from distutils.core import setup
import os
import shutil
shutil.copy('README.md', 'twochanradialwell/README.md')

dir_setup = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(dir_setup, 'twochanradialwell', 'release.py')) as f:
    # Defines __version__
    exec(f.read())

setup(name='twochanradialwell',
      version=__version__,
      description='Calculates solutions to the two channel radial well as described in Newton\'s "Scattering Theory of Waves and Particles".',
      author="Peter Bingham",
      author_email="petersbingham@hotmail.co.uk",
      packages=['twochanradialwell'],
      package_data={'twochanradialwell': ['tests/*', 'README.md']}
     )
