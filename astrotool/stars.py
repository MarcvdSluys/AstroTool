#!/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: EUPL-1.2

#  Copyright (c) 2019-2025  Marc van der Sluys - marc.vandersluys.nl
#  
#  This file is part of the AstroTool Python package,
#  see: https://www.nikhef.nl/~sluys/AstroTool/
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
    
    mass = _np.asarray(_np.copy(mass), dtype=float)    # Copy and typecast to numpy.ndarray
    scalar_input = False
    if mass.ndim == 0:
        mass = mass[None]  # Makes x a 1D array.  Comment: use np.newaxis instead?
        scalar_input = True
    
    rad = mass                                  # R/Ro ~ M/Mo for M <= Mo
    rad[mass>1] = _np.power(mass[mass>1], 0.6)  # R/Ro ~ (M/Mo)^0.6 for M >= Mo
    
    if scalar_input:
        return float(_np.squeeze(rad))  # Array -> scalar (float)
    
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


def rgb_coremass_at_R25Ro_from_mass(mass):
    """Return the helium core mass where an RGB star reaches a radius of 25Ro, as a function of the total
    mass.
    
    Parameters:
      mass (float):  Total mass (0.8-3.0Mo).
    
    Returns:
      (float):  Helium core mass.
    
    Notes:
      - Fit made for 0.8-3.0Mo, Z=0.02.
    """
    
    mass = _np.asarray(_np.copy(mass), dtype=float)  # Copy and typecast to numpy.ndarray
    scalar_input = False
    if mass.ndim == 0:
        mass = mass[None]  # Makes mass a 1D array.  Comment: use np.newaxis instead?
        scalar_input = True
    
    mc = _np.zeros_like(mass, dtype=float)  # Ensure float, calculation below goes wrong as int!
    mc[mass<=2.05] = 0.295 + 1.69e-3 * _np.power(mass[mass<=2.05], 4.51)  # Mean/med. |abs. diff.|: 0.000666 / 0.000562 Mo;  mean/med. |rel. diff.|: 0.00208 / 0.00159 (frac)
    mc[mass >2.05] = 0.278 + 2.13e-3 * _np.power(mass[mass >2.05], 3.79)  # Mean/med. |abs. diff.|: 0.000736 / 0.000452 Mo;  mean/med. |rel. diff.|: 0.00230 / 0.00139 (frac)
    
    if scalar_input:
        return float(_np.squeeze(mc))  # Arrays -> "scalars".  Note: type will still be np.array(scalar)
    
    return mc


def rgb_radius_from_mass_coremass(mass, mc):
    """Return an estimate of the radius of a low-mass (0.8-3Mo) RGB star from its total mass and He core mass.
    
    Parameters:
      mass (float):  Total stellar mass (Mo).
      mc (float):    Helium core mass (Mo).
    
    Returns:
      (float):  Radius of the star (Ro).
    
    Notes:
      - Fit made for 0.8-3.0Mo, Z=0.02.  Accuracy:
        - 0.8-1.3Mo: ~14-17%;
        - 1.4-2.0Mo: ~21-29%;
        - 2.1-2.3Mo: ~18-21%;
        - 2.4-2.6Mo: ~35-39%;
        - 2.7-3.0Mo: ~50-66%.
        - sticking to 0.8-2.0/2.3Mo can push the accuracy towards 10%:
          - Mmax = 3.0Mo: mean accuracy = 19% (more data points on longer low-mass tracks)
          - Mmax = 2.5Mo: mean accuracy = 11%
          - Mmax = 2.0Mo: mean accuracy =  8%
          - for more massive stars (>2.0-2.5Mo) a different fit function may be preferred
    """
    
    mcl = mc + rgb_coremass_at_R25Ro_from_mass(1) - rgb_coremass_at_R25Ro_from_mass(mass)  # Shift in core mass
    rad = 1.71961 + 1.96840e-3 * _np.exp(43.1952*mcl - 40.6587*_np.square(mcl))  # a + b*exp[c*x + d*x^2]
    return rad


def wd_radius_from_mass(mass):
    """Return the radius of a white dwarf from its mass.
    
    Parameters:
      mass (float):  White-dwarf mass (Mo).
    
    Returns:
      (float):  White-dwarf radius (Ro).
    
    Notes:
      - Nauenberg (1972), Eq.27 (https://ui.adsabs.harvard.edu/abs/1972ApJ...175..417N).
    """
    
    rad = _np.power(mass/1.454, 2/3)
    radius = 0.01125*_np.sqrt(1/rad - rad)  # Orig Nauenberg (1972), Eq.27
    
    # rad = np.power(df.mass/1.44, 2/3)
    # radius = 1.14e-2*np.sqrt(1/rad - rad)  # FV
    
    return radius


# Test code:
if __name__ == '__main__':
    pass
    
