# -*- coding: utf-8 -*-

from distutils.core import setup
import shutil
shutil.copy('README.md', 'twochanradialwell/README.md')

setup(name='twochanradialwell',
      version='0.13',
      description='Calculates solutions to the two channel radial well as described in Newton\'s "Scattering Theory of Waves and Particles".',
      author="Peter Bingham",
      author_email="petersbingham@hotmail.co.uk",
      packages=['twochanradialwell'],
      package_data={'twochanradialwell': ['tests/*', 'README.md']}
     )
