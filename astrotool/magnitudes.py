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


"""Magnitude-flux functions for AstroTool."""


# Allow relative imports from __main__() when running this file (PEP 366):
if __name__ == '__main__' and __package__ is None:
    __package__ = 'astrotool'

# Modules:
import numpy as _np

# Constants for each photometric band:
_band_mag0 = {'bol':18.99,                                            # Bolometric
              'U':25.90, 'B':25.36, 'V':26.02, 'R':26.66, 'I':27.37}  # Johnson UBVRI - Verbunt: Het leven van sterren, appendix
# 'GG':25.6884, 'GBP':25.3514, 'GRP':24.7619}             # Gaia: https://dc.g-vo.org/tableinfo/gaia.dr2epochflux
               
_band_eff_widths = {'bol':1,                                             # Bolometric: -
                    'U':68, 'B':98, 'V':89, 'R':138, 'I':149}            # Johnson UBVRI: Δλ in nm - Verbunt: Het leven van sterren, appendix


def flux_from_magnitude(mag, band='V'):
    """Convert magnitude to flux.
    
    Parameters:
      mag (float):   Magnitude of the object (-; ~ -2.5 log f)
      band (float):  Photometric band.
    
    Available photometric bands:
      - 'bol':  bolometric
      - 'U':    Johnson U
      - 'B':    Johnson B
      - 'V':    Johnson V
      - 'R':    Johnson R
      - 'I':    Johnson I
      - 'GG':   Gaia G   - not implemented!
      - 'GBP':  Gaia BP  - not implemented!
      - 'GRP':  Gaia RP  - not implemented!
    
    Returns:
      (float):  Flux (W/m2).
    """
    
    if band not in _band_mag0:
        print('flux_from_magnitude(): unknown band: ', band)
        print('Allowed bands are: ', _band_mag0.keys())
        print('Aborting...')
        exit(1)
        
    flux = _np.power(10, -(mag+_band_mag0[band])/2.5) * _band_eff_widths[band]
    
    return flux


def magnitude_from_flux(flux, band='V'):
    """Convert flux to magnitude.
    
    Parameters:
      flux (float):  Flux of the object (W/m2)
      band (float):  Photometric band.
    
    Available photometric bands:
      - 'bol':  bolometric
      - 'U':    Johnson U
      - 'B':    Johnson B
      - 'V':    Johnson V
      - 'R':    Johnson R
      - 'I':    Johnson I
      - 'GG':   Gaia G   - not implemented!
      - 'GBP':  Gaia BP  - not implemented!
      - 'GRP':  Gaia RP  - not implemented!
    
    Returns:
      (float):  Magnitude of the object (-; ~ -2.5 log f).
    """
    
    if band not in _band_mag0:
        print('magnitude_from_flux(): unknown band: ', band)
        print('Allowed bands are: ', _band_mag0.keys())
        print('Aborting...')
        exit(1)
        
    mag = -2.5 * _np.log10(flux) / _band_eff_widths[band] - _band_mag0[band]
    
    return mag


# Test code:
if __name__ == '__main__':
    import colored_traceback as _clrtrb
    _clrtrb.add_hook()
    
    print('Flux Mb=0  (2.54e-8  watt/m2/nm):  %10.3e' % (flux_from_magnitude(0.0, 'bol')))
    
    print('Flux U=0   (4.35e-11 watt/m2/nm):  %10.3e' % (flux_from_magnitude(0.0, 'U')))
    print('Flux B=0   (7.19e-11 watt/m2/nm):  %10.3e' % (flux_from_magnitude(0.0, 'B')))
    print('Flux V=0   (3.92e-11 watt/m2/nm):  %10.3e' % (flux_from_magnitude(0.0, 'V')))
    print('Flux R=0   (2.18e-11 watt/m2/nm):  %10.3e' % (flux_from_magnitude(0.0, 'R')))
    print('Flux I=0   (1.13e-11 watt/m2/nm):  %10.3e' % (flux_from_magnitude(0.0, 'I')))
    
    # print()
    # print('Flux Gaia G=0:   %10.3e' % (flux_from_magnitude(0.0, 'GG')))
    # print('Flux Gaia BP=0:  %10.3e' % (flux_from_magnitude(0.0, 'GBP')))
    # print('Flux Gaia RP=0:  %10.3e' % (flux_from_magnitude(0.0, 'GRP')))
    
    
    print()
    print()
    print('Flux -> Mb:  %6.3f' % (magnitude_from_flux(2.535e-8,  'bol')))
    
    print('Flux -> U:   %6.3f' % (magnitude_from_flux(4.365e-11, 'U')))
    print('Flux -> B:   %6.3f' % (magnitude_from_flux(7.178e-11, 'B')))
    print('Flux -> V:   %6.3f' % (magnitude_from_flux(3.908e-11, 'V')))
    print('Flux -> R:   %6.3f' % (magnitude_from_flux(2.168e-11, 'R')))
    print('Flux -> I:   %6.3f' % (magnitude_from_flux(1.127e-11, 'I')))
    
    # print()
    # print('Flux -> Gaia G:   %6.3f' % (magnitude_from_flux(5.304e-11, 'GG')))
    # print('Flux -> Gaia BP:  %6.3f' % (magnitude_from_flux(7.235e-11, 'GBP')))
    # print('Flux -> Gaia RP:  %6.3f' % (magnitude_from_flux(1.245e-10, 'GRP')))
    
    # print()
    # print('Flux I=0   (1.13e-11 watt/m2/nm):  %10.2e' % (magnitude_from_flux(0.0, 'BLA')))
 
