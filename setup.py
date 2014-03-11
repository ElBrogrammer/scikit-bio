#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, bipy development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

__version__ = '0.0.0-dev'

from setuptools import find_packages, setup
from distutils.command.build_py import build_py

classes = """
    Development Status :: 1 - Planning
    License :: OSI Approved :: BSD License
    Topic :: Software Development :: Libraries
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python
    Programming Language :: Python :: 2.7
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

long_description = """The bipy project"""

setup(name='bipy',
      cmdclass={'build_py': build_py},
      version=__version__,
      license='BSD',
      description='bipy',
      long_description=long_description,
      author="bipy development team",
      author_email="gregcaporaso@gmail.com",
      maintainer="bipy development team",
      maintainer_email="gregcaporaso@gmail.com",
      url='https://github.com/biocore/bipy',
      packages=find_packages(),
      install_requires=['numpy >= 1.5.1', 'matplotlib >= 1.1.0',
                        'scipy >=0.13.0', 'sphinx-bootstrap-theme', 'numpydoc'],
      extras_require={'test': ["nose >= 0.10.1", "pep8"],
                      'doc': ["Sphinx >= 1.1", "numpydoc",
                              "sphinx-bootstrap-theme"]},
      classifiers=classifiers
      )
