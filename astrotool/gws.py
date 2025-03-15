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


"""Functions to deal with gravitational waves for AstroTool."""


# Allow relative imports from __main__() when running this file (PEP 366):
if __name__ == '__main__' and __package__ is None:
    __package__ = 'astrotool'

# Modules:
import numpy as _np
import astroconst as _ac


def chirp_mass_from_m1_m2(M1,M2):
    """Compute the chirp mass from the individual masses.
    
    Parameters:
      M1 (float):  Mass 1.
      M2 (float):  Mass 2 (same unit as M1).
    
    Returns:
      (float):  Chirp mass (same unit as M1 and M2).
    """
    
    Mtot3  = _np.power(M1+M2, 1/3)
    Mchirp = _np.power(M1 * M2 / Mtot3, 3/5)
    
    return Mchirp


def time_in_frequency_band(Mc, fmin,fmax):
    """Compute the time a GW source/signal with given chirp mass spends in a given frequency band.
    
    Parameters:
      Mc (float):    Chirp mass of the source (Mo).
      fmin (float):  Lower limit of the frequency band (Hz).
      fmax (float):  Upper limit of the frequency band (Hz).
    
    Returns:
      (float):  Time the source spends in the frequency band (s).
    """
    
    dt_obs = 5/256 * _ac.pi**(-8/3) * _np.power(_ac.G * Mc * _ac.Mo / _ac.c**3, -5/3) * \
        (fmin**(-8/3) - fmax**(-8/3))
    
    return dt_obs


def strain_optimal_amplitude_from_masses_and_orbit(M1,M2, aorb, dist):
    """Compute the optimal strain amplitude from the individual masses, orbital separation and distance.
    
    Parameters:
      M1 (float):    Mass 1 (Mo).
      M2 (float):    Mass 2 (Mo).
      aorb (float):  Orbital separation (Ro).
      dist (float):  Distance (Mpc).
    
    Returns:
      (float):  Optimal strain amplitude h_+ or h_x (-).
    """
    
    h_opt = 4/(dist * _ac.Mpc)  *  _ac.G**2 / _ac.c**4  *  M1 * M2 * _ac.Mo**2  /  (aorb * _ac.Ro)
    
    return h_opt

    
def strain_optimal_amplitude_from_masses_and_frequency(M1,M2, f_gw, dist):
    """Compute the optimal strain amplitude from the individual masses, GW frequency and distance.
    
    Parameters:
      M1 (float):    Mass 1 (Mo).
      M2 (float):    Mass 2 (Mo).
      f_gw (float):  GW frequency (Hz) (=2x f_orb!).
      dist (float):  Distance (Mpc).
    
    Returns:
      (float):  Optimal strain amplitude h_+ or h_x (-).
    """
    
    Mtot3 = _np.power((M1+M2)*_ac.Mo, 1/3)
    h_opt = 4/(dist * _ac.Mpc)  *  _ac.G**(5/3) / _ac.c**4  *  M1 * M2 * _ac.Mo**2  /  Mtot3 \
        * (_ac.pi * f_gw)**(2/3)
    
    return h_opt

    
def strain_optimal_amplitude_from_chirpmass_and_frequency(Mc, f_gw, dist):
    """Compute the optimal strain amplitude from the chirp mass, GW frequency and distance.
    
    Parameters:
      Mc (float):    Chirp mass (Mo).
      f_gw (float):  GW frequency (Hz) (=2x f_orb!).
      dist (float):  Distance (Mpc).
    
    Returns:
      (float):  Optimal strain amplitude h_+ or h_x (-).
    """
    
    h_opt = 4/(dist * _ac.Mpc)  *  _np.power(_ac.G * Mc * _ac.Mo / _ac.c**2, 5/3)  \
        * (_ac.pi * f_gw / _ac.c)**(2/3)
    
    return h_opt

    
def time_to_merger_from_a(m1,m2, aorb):
    """Return the time to merger due to GWs for a binary with given masses and orbital separation.
    
    Parameters:
      m1 (float):    Mass of star 1 (Mo).
      m2 (float):    Mass of star 2 (Mo).
      aorb (float):  Orbital separation (Ro).
    
    Returns:
      (float):  Merger time (yr).
    """
    from .binaries import gw_merger_time_from_a as _gw_merger_time_from_a
    return _gw_merger_time_from_a(m1,m2, aorb)
    

def time_to_merger_from_P(m1,m2, Porb):
    """Return the merger time due to GWs for a binary with given masses and orbital period.
    
    Parameters:
      m1 (float):    Mass of star 1 (Mo).
      m2 (float):    Mass of star 2 (Mo).
      Porb (float):  Orbital period (days).
    
    Returns:
      (float):  Merger time (yr).
    """
    from .binaries import gw_merger_time_from_P as _gw_merger_time_from_P
    return _gw_merger_time_from_P(m1,m2, Porb)
    

def isco_gw_frequency_from_mass(mass):
    """Return the ISCO GW frequency for a black hole with given mass.
    
    Parameters:
      mass (float):  Black-hole mass (Mo).
    
    Returns:
      (float):  GW ISCO frequency (Hz)  (=2x f_orb!).
    """
    
    return 1/(6**(3/2)*_ac.pi) * _ac.c**3 / (_ac.G * mass * _ac.Mo)


# Test code:
if __name__ == '__main__':
    pass
    
