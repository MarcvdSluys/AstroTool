#!/bin/env python3
# -*- coding: utf-8 -*-

"""Setup.py for the SolarEnergy Python package."""


# Package version:
version='0.0.9'

# Get long description from README.md:
with open('README.md', 'r') as fh:
    long_description = fh.read()


from setuptools import setup
setup(
    name='astrotool',
    description='A Python package for astronomical calculations in Python or on the command line',
    author='Marc van der Sluys',
    url='http://astro.ru.nl/~sluys/AstroTool',
    
    packages=['astrotool'],
    install_requires=['astroconst','colored_traceback','numpy'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    
    version=version,
    license='EUPL 1.2',
    keywords=['astronomy','ephemeris','date','time','coordinates'],
    
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: European Union Public Licence 1.2 (EUPL 1.2)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Astronomy',
    ]
)
