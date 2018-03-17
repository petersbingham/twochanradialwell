# -*- coding: utf-8 -*-

from distutils.core import setup
import shutil
shutil.copy('README.md', 'TwoChanRadialWell/README.md')

setup(name='TwoChanRadialWell',
      version='0.5',
      description='Calculates solutions to the two channel radial well as described in Newton\'s "Scattering Theory of Waves and Particles".',
      author="Peter Bingham",
      author_email="petersbingham@hotmail.co.uk",
      packages=['TwoChanRadialWell'],
      package_data={'TwoChanRadialWell': ['tests/*', 'README.md']}
     )
