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
import astroconst as _ac

# Constants for each photometric band: f_0 in W/m^2:
_band_flux0 = {'bol':2.518e-8,                                                   # Bolometric - IAU Resolution B2 (2015): https://arxiv.org/abs/1510.06262
               'U':3.97e-9, 'B':6.13e-9, 'V':3.63e-9, 'R':2.17e-9, 'I':1.13e-9,  # Bessell/Johnson UBVRI: https://svo2.cab.inta-csic.es/theory/fps/index.php?gname=Generic&gname2=Bessell
               'Ks_SOFI':4.21e-11,                                               # Ks (IR) filter on La Silla/NTT/SOFI: https://svo2.cab.inta-csic.es/theory/fps/index.php?gname=LaSilla&gname2=SOFI
               'GG':2.49e-9, 'GBP':4.04e-9, 'GRP':1.28e-9,                       # Gaia DR2 revised filters: https://svo2.cab.inta-csic.es/theory/fps/index.php?gname=GAIA&gname2=GAIA2r
               'SDSS_u':3.75e-9,   'SDSS_g':5.45e-9,   'SDSS_r':2.5e-9,    'SDSS_i':1.39e-9,   'SDSS_z':8.39e-10,   # SDSS filters: https://svo2.cab.inta-csic.es/theory/fps/index.php?gname=SLOAN
               'SDSS_upr':3.56e-9, 'SDSS_gpr':5.28e-9, 'SDSS_rpr':2.42e-9, 'SDSS_ipr':2.42e-9, 'SDSS_zpr':4.88e-10  # SDSS prime filters: https://svo2.cab.inta-csic.es/theory/fps/index.php?gname=SLOAN
               }


def flux_from_magnitude(mag, band='V'):
    """Convert magnitude to flux.
    
    Parameters:
      mag (float):   Magnitude of the object (-; ~ -2.5 log f)
      band (float):  Photometric band.
    
    Available photometric bands:
      - 'bol':                Bolometric
      - 'U','B','V','R','I':  Johnson UBVRI
      - 'Ks_SOFI:             Ks (IR) filter on La Silla/NTT/SOFI
      - 'GG','GBP','GRP':     Gaia G, BP, RP
    
    Returns:
      (float):  Flux (W/m2).
    """
    
    if band not in _band_flux0:
        print('flux_from_magnitude(): unknown band: ', band)
        print('Allowed bands are: ', _band_flux0.keys())
        print('Aborting...')
        exit(1)
        
    flux = _band_flux0[band] * _np.power(10, -mag/2.5)
    
    return flux


def magnitude_from_flux(flux, band='V'):
    """Convert flux to magnitude.
    
    Parameters:
      flux (float):  Flux of the object (W/m2)
      band (float):  Photometric band.
    
    Available photometric bands:
      - 'bol':                Bolometric
      - 'U','B','V','R','I':  Johnson UBVRI
      - 'Ks_SOFI:             Ks (IR) filter on La Silla/NTT/SOFI
      - 'GG','GBP','GRP':     Gaia G, BP, RP
    
    Returns:
      (float):  Magnitude of the object (-; ~ -2.5 log f).
    """
    
    if band not in _band_flux0:
        print('magnitude_from_flux(): unknown band: ', band)
        print('Allowed bands are: ', _band_flux0.keys())
        print('Aborting...')
        exit(1)
        
    mag = -2.5 * _np.log10(flux/_band_flux0[band])
    
    return mag


def luminosity_from_flux_at_10pc(flux):
    """Compute the luminosity in W from the flux measured at a distance of 10pc.
    
    Useful after conversion of absolute magnitude to flux.
    
    Parameters:
      flux (float):  Flux of the object in some band or bolometric (W/m2).
    
    Returns:
      (float):  Luminosity of the object in that band (W).
    """
    
    lum = flux * 4*_ac.pi*(10*_ac.pc)**2
    
    return lum


# Test code:
if __name__ == '__main__':
    import colored_traceback as _clrtrb
    _clrtrb.add_hook()
    
    print('Flux Mb=0  (2.518e-8 watt/m2):  %10.3e' % (flux_from_magnitude(0.0, 'bol')))
    
    print('Flux U=0   (3.97e-9  watt/m2):  %10.3e' % (flux_from_magnitude(0.0, 'U')))
    print('Flux B=0   (6.13e-9  watt/m2):  %10.3e' % (flux_from_magnitude(0.0, 'B')))
    print('Flux V=0   (3.63e-9  watt/m2):  %10.3e' % (flux_from_magnitude(0.0, 'V')))
    print('Flux R=0   (2.17e-9  watt/m2):  %10.3e' % (flux_from_magnitude(0.0, 'R')))
    print('Flux I=0   (1.13e-9  watt/m2):  %10.3e' % (flux_from_magnitude(0.0, 'I')))
    
    print()
    print('Flux Ks=0  (4.21e-11 watt/m2):  %10.3e' % (flux_from_magnitude(0.0, 'Ks_SOFI')))
    
    print()
    print('Flux Gaia G=0:   %10.3e' % (flux_from_magnitude(0.0, 'GG')))
    print('Flux Gaia BP=0:  %10.3e' % (flux_from_magnitude(0.0, 'GBP')))
    print('Flux Gaia RP=0:  %10.3e' % (flux_from_magnitude(0.0, 'GRP')))
    
    
    print()
    print()
    print('Flux -> Mb:  %6.3f' % (magnitude_from_flux(2.518e-8,  'bol')))
    
    print('Flux -> U:   %6.3f' % (magnitude_from_flux(3.97e-9, 'U')))
    print('Flux -> B:   %6.3f' % (magnitude_from_flux(6.13e-9, 'B')))
    print('Flux -> V:   %6.3f' % (magnitude_from_flux(3.63e-9, 'V')))
    print('Flux -> R:   %6.3f' % (magnitude_from_flux(2.17e-9, 'R')))
    print('Flux -> I:   %6.3f' % (magnitude_from_flux(1.13e-9, 'I')))
    
    print()
    print('Flux -> SOFI Ks:  %6.3f' % (magnitude_from_flux(4.21e-11, 'Ks_SOFI')))
    
    print()
    print('Flux -> Gaia G:   %6.3f' % (magnitude_from_flux(2.49e-9, 'GG')))
    print('Flux -> Gaia BP:  %6.3f' % (magnitude_from_flux(4.04e-9, 'GBP')))
    print('Flux -> Gaia RP:  %6.3f' % (magnitude_from_flux(1.28e-9, 'GRP')))
    
    # print()
    # print('Flux I=0   (1.13e-11 watt/m2/nm):  %10.2e' % (magnitude_from_flux(0.0, 'BLA')))
 
    print()
    print('Flux @10pc -> L(W):  %0.5g' % (luminosity_from_flux_at_10pc(1e-10)))  # 1.1965e+26
    
