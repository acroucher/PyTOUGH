# PyTOUGH setup script
from distutils.core import setup

setup(name = 'PyTOUGH',
      version = '1.5.1',
      description = 'Python scripting library for TOUGH2 simulation',
      author = 'Adrian Croucher',
      author_email = 'a.croucher@auckland.ac.nz',
      url = 'https://github.com/acroucher/PyTOUGH',
      license = 'LGPL',
      py_modules = [
          'fixed_format_file',
          'geometry',
          'IAPWS97',
          'mulgrids',
          't2data',
          't2grids',
          't2incons',
          't2listing',
          't2thermo'],
      )
