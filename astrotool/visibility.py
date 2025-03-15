#!/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: EUPL-1.2

#  Copyright (c) 2019-2025  Marc van der Sluys - marc.vandersluys.nl
#   
#  This file is part of the AstroTool Python package,
#  see: https://www.nikhef.nl/~sluys/AstroTool/
#   
#  AstroTool has been developed by Marc van der Sluys of the Department of Physics at Utrecht
#  University in the Netherlands, and the Netherlands Institute for Nuclear and High-Energy Physics (Nikhef)
#  in Amsterdam.
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


"""Object visibility functions for AstroTool."""


# Allow relative imports from __main__() when running this file (PEP 366):
if __name__ == '__main__' and __package__ is None:
    __package__ = 'astrotool'

# Modules:
import numpy as _np
import astroconst as _ac


def sky_brightness_cdm2_from_nL(nL):
    """Convert sky brightness in nanolambert to candela per square meter.
    
    Parameters:
      nL (float):  Sky brightness in nanolambert (nL).
    
    Returns:
      (float):  Sky brightness in candela per square meter (cd/m^2).
    """
    
    return nL * 1e-5/_ac.pi


def sky_brightness_nL_from_cdm2(cdm2):
    """Convert sky brightness in candela per square meter to nanolambert.
    
    Parameters:
      cdm2 (float):  Sky brightness in candela per square meter (cd/m^2).
    
    Returns:
      (float):  Sky brightness in nanolambert (nL).
    """
    
    return cdm2 * 1e5 * _ac.pi


def sky_brightness_cdm2_from_mas2(mas2):
    """Convert sky brightness in magnitude per square arcsecond to candela per square meter.
    
    Parameters:
      mas2 (float):  Sky brightness in magnitude per square arcsecond.
    
    Returns:
      (float):  Sky brightness in candela per square meter.
    
    Note: conversion seems to come from the Sky Quality Meter (SQM).
    """
    
    return 1.08e5 * 10**(-mas2/2.5)  # mag/arcsec^2 -> cd/m^2


def sky_brightness_mas2_from_cdm2(cdm2):
    """Convert sky brightness in candela per square meter to magnitude per square arcsecond.
    
    Parameters:
      cdm2 (float):  Sky brightness in candela per square meter.
    
    Returns:
      (float):  Sky brightness in magnitude per square arcsecond.
    
    Note: conversion seems to come from the Sky Quality Meter (SQM).
    """
    
    return -2.5*_np.log10(cdm2/1.08e5)  # mag/arcsec^2 -> cd/m^2


def sky_brightness_mlim_from_cdm2(cdm2):
    """Convert sky brightness in candela per square meter to a limiting visual magnitude.
    
    Parameters:
      cdm2 (float):  Sky brightness in cd/m^2.
    
    Returns:
      (float):  Sky brightness as limiting magnitude.
    
    Note:
      Adapted from Weaver (1947) and Garsteng (1986):
      - shift to match data by Weaver (1947): 7.930 -> 7.706;  4.305 -> 4.115;
      - use candela per square meter instead of nanolambert.
    """
    
    cdm2 = _np.asarray(_np.copy(cdm2))     # Copy and typecast to numpy.ndarray
    if cdm2.ndim == 0:  cdm2 = cdm2[None]  # Makes cdm2 1D.
    
    nL = sky_brightness_nL_from_cdm2(cdm2)  # Convert from candela per square meter to nanolambert
    
    Mlim = _np.zeros_like(cdm2)  # Need an array of the same size
    Mlim[nL<=1479] = 7.706 - 5*_np.log10(1 + 0.1122   * _np.sqrt(nL[nL<=1479]))   # b <= 1479 nL = Mlim>4.0
    Mlim[nL>1479]  = 4.115 - 5*_np.log10(1 + 0.001122 * _np.sqrt(nL[nL>1479]))    # b >  1479 nL = Mlim<4.0
    
    return _np.squeeze(Mlim)


def sky_brightness_cdm2_from_mlim(mlim):
    """Convert limiting visual magnitude to sky brightness in candela per square meter.
    
    Parameters:
      mlim (float):  Sky brightness as limiting magnitude.
    
    Returns:
      (float):  Sky brightness in candela per square meter.
    
    Note:
      Inverse of sky_brightness_mlim_from_cdm2() - see documentation there.
    """
    
    mlim = _np.asarray(_np.copy(mlim))     # Copy and typecast to numpy.ndarray
    if mlim.ndim == 0:  mlim = mlim[None]  # Makes mlim 1D.
    
    nL = _np.zeros_like(mlim)  # Need an array of the same size - nanolambert
    nL[mlim>=4] = _np.square((10**(-(mlim[mlim>=4] - 7.706)/5) - 1)/0.1122)    # Mlim >= 4.0; b <= 1479 nL
    nL[mlim<4]  = _np.square((10**(-(mlim[mlim<4]  - 4.115)/5) - 1)/0.001122)  # Mlim < 4.0;  b >  1479 nL
    
    cdm2 = sky_brightness_cdm2_from_nL(nL)  # Convert nanolambert to candela per square meter
    
    return _np.squeeze(cdm2)


# Test code:
if __name__ == '__main__':
    import colored_traceback as _clrtrb
    _clrtrb.add_hook()
    
