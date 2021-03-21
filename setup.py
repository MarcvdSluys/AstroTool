#!/bin/env python3

"""Setup.py for the SolarEnergy Python package."""


# Package version:
version="0.0.1"

# Get long description from README.md:
with open("README.md", "r") as fh:
    long_description = fh.read()


from setuptools import setup
setup(
    name='astrotool',
    description='A Python package for astronomical calculations in Python or on the command line',
    author='Marc van der Sluys',
    url='http://astro.ru.nl/~sluys/AstroTool',
    
    packages=['astrotool'],
    install_requires=['numpy','fortranformat'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    
    version=version,
    license='GPLv3+',
    keywords=['astronomy','ephemeris','date','time','coordinates'],
    
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Astronomy",
    ]
)
