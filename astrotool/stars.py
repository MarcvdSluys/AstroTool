#!/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: EUPL-1.2

#  Copyright (c) 2019-2024  Marc van der Sluys - marc.vandersluys.nl
#  
#  This file is part of the AstroTool Python package,
#  see: http://astro.ru.nl/~sluys/AstroTool/
#  
#  AstroTool has been developed by Marc van der Sluys of the Department of Physics at Utrecht
#  University in the Netherlands, and the Netherlands Institute for Nuclear and High-Energy
#  Physics (Nikhef) in Amsterdam.
#  
#  This is free software: you can redistribute it and/or modify it under the terms of the
#  European Union Public Licence 1.2 (EUPL 1.2).
#  
#  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
#  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the EU Public Licence for more details.
#  
#  You should have received a copy of the European Union Public Licence along with this code.
#  If not, see <https://www.eupl.eu/1.2/en/>.


"""Functions on stars and stellar evolution for AstroTool."""


# Allow relative imports from __main__() when running this file (PEP 366):
if __name__ == '__main__' and __package__ is None:
    __package__ = 'astrotool'

# Modules:
import numpy as _np
# import astroconst as _ac


def ms_rad_from_mass(mass):
    """Return an estimate of the radius of a main-sequence star from its mass (e.g. Verbunt 2010, Eq.4.5).
    
    Parameters:
      mass (float):  mass of the star (Mo).
    
    Returns:
      (float):  Radius of the star (Ro).
    
    Note: fit on lower bound of R (i.e. ~R_zams), using 183 stars from 96 double-lined eclipsing binaries and
    the Sun, for which M and R are known to an accuracy of 2% or better.  Valid for 0.12 ~< M/Mo ~< 25, not so
    great for M~2.5Mo.
    """
    
    mass = _np.asarray(_np.copy(mass))    # Copy and typecast to numpy.ndarray
    scalar_input = False
    if mass.ndim == 0:
        mass = mass[None]  # Makes x a 1D array.  Comment: use np.newaxis instead?
        scalar_input = True
    
    rad = mass                                  # R/Ro ~ M/Mo for M <= Mo
    rad[mass>1] = _np.power(mass[mass>1], 0.6)  # R/Ro ~ (M/Mo)^0.6 for M >= Mo
    
    if scalar_input:
        return _np.squeeze(rad)  # Arrays -> "scalars".  Note: type will still be _np.array(scalar)
    
    return rad


def ms_lum_from_mass(mass):
    """Return an estimate of the luminosity of a main-sequence star from its mass (e.g. Verbunt 2010, Eq.4.6).
    
    Parameters:
      mass (float):  mass of the star (Mo).
    
    Returns:
      (float):  Luminosity of the star (Lo).
    
    Note: fit using 183 stars from 96 double-lined eclipsing binaries and the Sun, for which M and L are known
    to an accuracy of 2% or better.  Valid for 0.2 ~< M/Mo ~< 25; pretty tight fit (in logM-logL).
    """
    
    return _np.power(mass, 3.8)


def rgb_rad_from_mc_it84(mc):
    """Return an estimate of the radius of a (red) giant based on its (He) core mass, using Iben and Tutukov
    (1984).
    
    Parameters:
      mc (float):  (He) core mass (Mo).
    
    Returns:
      (float):  Radius (Ro).
    """
    
    return 10**3.5 * mc**4


# Test code:
if __name__ == '__main__':
    pass
    
