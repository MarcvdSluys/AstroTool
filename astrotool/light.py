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


"""Functions to do with light for AstroTool."""


# Allow relative imports from __main__() when running this file (PEP 366):
if __name__ == '__main__' and __package__ is None:
    __package__ = 'astrotool'


# Modules:
import numpy as _np
import astroconst as _ac


def planck_flux_from_lambda(Teff, lam, per_steradian=False):
    """Compute the Planck spectrum for a black body with given temperature and for a given (range of)
    wavelength(s).
    
    Parameters:
      Teff (float):          Effective temperature (K)
      lam (float):           Wavelength (nm)
      per_steradian (bool):  Return result per steradian of solid angle (defaults to False)
    
    Returns:
      (float): Flux/intensity/irradiation in Watt / m^2 / nm  (/ steradian)
    """
    
    exp = _np.exp(_ac.h_p * _ac.c / (lam*_ac.nm * _ac.k_b * Teff))
    B_lam = 2 * _ac.h_p * _ac.c**2 / _np.power(lam*_ac.nm,5) * 1 / (exp - 1) * _ac.nm
    if not per_steradian: B_lam *= _ac.pi  # "Integrate" over solid angle
    
    return B_lam


def planck_flux_from_frequency(Teff, freq, per_steradian=False):
    """Compute the Planck function/black-body spectrum from frequency and effective temperature.
    
    Parameters:
      Teff (float):          Effective temperature (K)
      freq (float):          EM frequency (Hz)
      per_steradian (bool):  Return result per steradian of solid angle (defaults to False)
    
    Returns:
      (float):  Flux/intensity/irradiation in Watt / m^2 / Hz  (/ steradian)
    """
    
    exp = _ac.h_p * freq / (_ac.k_b * Teff)
    B_nu = 2 * _ac.h_p * _np.power(freq,3) / _ac.c**2 / (_np.exp(exp) - 1)
    if not per_steradian: B_nu *= _ac.pi  # "Integrate" over solid angle
    
    return B_nu


def planck_peak_lambda(Teff):
    """Compute the wavelength where the Planck spectrum peaks for a given effective temperature.
    
    Parameters:
      Teff (float):          Effective temperature (K)
    
    Returns:
      (float): Wavelength where the Planck spectum peaks (lambda_max) in nm.
    """
    
    lambda_max = 2891941.5 / Teff  # nm K -> nm;  0.201 h*c/k = 2891941.5 nm K
    
    return lambda_max


def planck_peak_freq(Teff):
    """Compute the frequency where the Planck spectrum peaks for a given effective temperature.
    
    Parameters:
      Teff (float):          Effective temperature (K)
    
    Returns:
      (float): Frequency where the Planck spectum peaks (nu_max) in Hz.
    """
    
    nu_max = 5.8759266e+10 * Teff  # Hz/K -> Hz;  2.82 k/h 5.8759266e+10 Hz/K
    
    return nu_max


# Test code:
if __name__ == '__main__':
    _TeffSun = 5771  # K
    _lam_max = planck_peak_lambda(_TeffSun)
    _nu_max  = planck_peak_freq(_TeffSun)
    
    print('Solar lambda_max:               %0.4f nm' % (_lam_max))
    print('Solar peak flux at lambda_max:  %0.4e W/m2/nm' % (planck_flux_from_lambda(_TeffSun, _lam_max)))
    print('Solar nu_max:                   %0.4e Hz' % (_nu_max))
    print('Solar peak flux at nu_max:      %0.4e W/m2/Hz' % (planck_flux_from_frequency(_TeffSun, _nu_max)))


