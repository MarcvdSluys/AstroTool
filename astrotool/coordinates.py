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


"""Coordinate transformations and related functions for AstroTool."""


# Allow relative imports from __main__() when running this file (PEP 366):
if __name__ == '__main__' and __package__ is None:
    __package__ = 'astrotool'

# Modules:
import numpy as _np
import astroconst as _ac
from .date_time import tjc_from_jd


def obliquity(jd):
    """Compute the obliquity of the ecliptic in radians from the JD(E).
    
    Args:
      jd (float):   Julian day (days).
    
    Returns:
      float:  eps: Obliquity of the ecliptic (rad).
    
    References:
      - Seidelman 1992, Eq. 3.222-1.

    """
    
    tJC = tjc_from_jd(jd)  # Time in Julian centuries since J2000.0
    eps = 0.409092804 - 2.269655e-4*tJC - 2.86e-9*tJC**2 + 8.78967e-9*tJC**3  # Obliquity of the ecliptic (rad)
    
    return eps


def eq2ecl(ra,dec, eps):
    """Convert equatorial coordinates to ecliptical.

    Args:
      ra (float):    Right ascension (rad).
      dec (float):   Declination (rad).
      eps (float):   Obliquity of the ecliptic (rad).
    
    Returns:
      tuple (float,float):  tuple containing (lon, lat):
    
        - lon (float):   Ecliptical longitude (rad).
        - lat (float):   Ecliptical latitude (rad).

    """
    
    lon = _np.arctan2( _np.sin(ra)  * _np.cos(eps) + _np.tan(dec) * _np.sin(eps),  _np.cos(ra) ) % _ac.pi2
    lat =  _np.arcsin( _np.sin(dec) * _np.cos(eps) - _np.cos(dec) * _np.sin(eps) * _np.sin(ra) )
    
    return lon,lat


def ecl2eq(lon,lat, eps):
    """Convert (geocentric) spherical ecliptical coordinates to spherical equatorial coordinates.
    
    Args:
      lon (float):   Ecliptical longitude (rad).
      lat (float):   Ecliptical latitude (rad).
      eps (float):   Obliquity of the ecliptic (rad).
    
    Returns:
      tuple (float,float):  tuple containing (ra, dec):
    
        - ra (float):    Right ascension (rad).
        - dec (float):   Declination (rad).
    
    References:
      - `Explanatory Supplement to the Astronomical Almanac 3rd Ed,
        Eq.14.43 <https://aa.usno.navy.mil/publications/docs/exp_supp.php>`_

    """
    
    ra  = _np.arctan2( _np.sin(lon) * _np.cos(eps)  -  _np.tan(lat) * _np.sin(eps),  _np.cos(lon) ) % _ac.pi2
    dec =  _np.arcsin( _np.sin(lat) * _np.cos(eps)  +  _np.cos(lat) * _np.sin(eps) * _np.sin(lon) )
    
    return ra,dec


def par2horiz(ha,dec, phi):
    """Convert parallactic coordinates to horizontal.
    
    Args:
      ha (float):    Hour angle (rad).
      dec (float):   Declination (rad).
      phi (float):   Geographical latitude (rad, N>0).
    
    Returns:
      tuple (float,float):  tuple containing (az, alt):
    
        - az (float):    Azimuth (rad, S=0).
        - alt (float):   Altitude (rad).
    
    """
    
    az  = _np.arctan2( _np.sin(ha),   _np.cos(ha) * _np.sin(phi) - _np.tan(dec) * _np.cos(phi) ) % _ac.pi2
    alt = _np.arcsin(  _np.sin(dec) * _np.sin(phi) + _np.cos(ha) * _np.cos(dec) * _np.cos(phi) )
    
    return az,alt


def proper_motion(jd_start,jd_target, ra,dec, pma,pmd):
    """Compute the proper motion from jd_start to jd_target for the given positions and proper motions.
    
    Args:
      jd_start (float):    Julian day of the initial epoch (days).
      jd_target (float):   Julian day of the target epoch (days).
    
      ra (float):          Right ascension (numpy array, rad).
      dec (float):         Declination (numpy array, rad).
    
      pma (float):         Proper motion in right ascension (numpy array, rad/yr).
      pmd (float):         Proper motion in declination (numpy array, rad/yr).
    
    Returns:
      tuple (float,float):  tuple containing (ra_target, dec_target):
    
        - ra_target (float):    Right ascension for the target epoch (rad).
        - dec_target (float):   Declination for the target epoch (rad).
    
    """
    
    dtYr   = (jd_target - jd_start)/365.25
    ra_target  = ra  + pma*dtYr / _np.cos(dec)
    dec_target = dec + pmd*dtYr
    
    return ra_target,dec_target


def precess_from_2000(jd, ra,dec):
    """Compute precession in equatorial coordinates from J2000 to that of the specified JD.
    
    J2000 is the equinox of many catalogues, including the Hipparcos one.
    
    Args:
      jd (float):    Julian day (days).
      ra (float):    Right ascension (rad).
      dec (float):   Declination (rad).
    
    Returns:
      tuple (float,float):  tuple containing (ra_new, dec_new):
    
        - ra_new (float):    Right ascension for the target equinox (rad).
        - dec_new (float):   Declination for the target equinox (rad).
    
    """
    
    tJC  = tjc_from_jd(jd)  # Time in Julian centuries since J2000.0
    tJC2 = tJC**2
    tJC3 = tJC*tJC2
    
    zeta  = (2306.2181*tJC + 0.30188*tJC2 + 0.017998*tJC3)*_ac.as2r  # Convert the result from arcseconds to radians
    z     = (2306.2181*tJC + 1.09468*tJC2 + 0.018203*tJC3)*_ac.as2r
    theta = (2004.3109*tJC - 0.42665*tJC2 - 0.041833*tJC3)*_ac.as2r
    
    ra_new  = (_np.arctan2( _np.sin(ra + zeta) * _np.cos(dec),  _np.cos(ra + zeta) * _np.cos(theta) * _np.cos(dec) - _np.sin(theta) * _np.sin(dec) ) + z) % _ac.pi2
    dec_new = _np.arcsin( _np.cos(ra + zeta) * _np.sin(theta) * _np.cos(dec)  +  _np.cos(theta) * _np.sin(dec) )
    
    return ra_new,dec_new


def geoc2topoc_ecl(lon_gc,lat_gc, dist_gc,rad_gc, eps,lst, lat_obs,ele_obs=0, debug=False):
    """Convert spherical ecliptical coordinates from the geocentric to the topocentric system.
    
    Args:
      lon_gc (float):    Geocentric ecliptic longitude (rad).
      lat_gc (float):    Geocentric ecliptic latitude (rad).
      dist_gc (float):   Geocentric distance (AU).
      rad_gc (float):    Geocentric semi-diameter (rad).
      
      eps (float):       Obliquity of the ecliptic (rad).
      lst (float):       Local sidereal time (rad).
      
      lat_obs (float):   Geographical latitude of the observer (rad).
      ele_obs (float):   Altitude/elevation of the observer above sea level (metres, optional, default value = 0).
      
      debug (bool):      Print debug output (True/False, optional, default value = True).
    
    Returns:
      tuple (float,float,float):  tuple containing (lon_tc, lat_tc, rad_tc):
    
        - lon_tc (float):   Topocentric ecliptic longitude (rad).
        - lat_tc (float):   Topocentric ecliptic latitude (rad).
        - rad_tc (float):   Topocentric semi-diameter (rad).
    
    """
    
    # Meeus, Ch.11, p.82:
    Req = _ac.r_earth        # Equatorial radius of the Earth in metres (same units as the elevation)
    #                          (http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf)
    RpolEq = 0.996647189335  # Rpol/Req = 1-f: flattening of the Earth - WGS84 ellipsoid 
    
    u  = _np.arctan(RpolEq*_np.tan(lat_obs))
    RsinPhi = RpolEq*_np.sin(u) + ele_obs/Req * _np.sin(lat_obs)
    RcosPhi = _np.cos(u)        + ele_obs/Req * _np.cos(lat_obs)
    
    sinHp = _np.sin(_ac.r_earth/_ac.au) / (dist_gc/_ac.au)  # Sine of the horizontal parallax, Meeus, Eq. 40.1
    
    # Meeus, Ch.40, p.282:
    N  = _np.cos(lon_gc)*_np.cos(lat_gc) - RcosPhi*sinHp*_np.cos(lst)
    
    lon_tc = _np.arctan2( _np.sin(lon_gc)*_np.cos(lat_gc) - sinHp*(RsinPhi*_np.sin(eps) + RcosPhi*_np.cos(eps)*_np.sin(lst)), N ) % _ac.pi2  # Topocentric longitude
    lat_tc = _np.arctan((_np.cos(lon_tc)*(_np.sin(lat_gc) - sinHp*(RsinPhi*_np.cos(eps) - RcosPhi*_np.sin(eps)*_np.sin(lst))))/N)         # Topocentric latitude
    rad_tc = _np.arcsin(_np.cos(lon_tc)*_np.cos(lat_tc)*_np.sin(rad_gc)/N)                                                                # Topocentric semi-diameter
    
    # print(dist_gc, dist_gc*rad_gc/rad_tc)
    
    if debug:
        print()
        print('geoc2topoc_ecl():')
        
        print('\nInput:')
        print('lon_gc: ', lon_gc)
        print('lat_gc: ', lat_gc)
        print('dist_gc: ', dist_gc)
        print('rad_gc: ', rad_gc)
        print('eps: ', eps)
        print('lst: ', lst)
        print('lat_obs: ', lat_obs)
        print('ele_obs: ', ele_obs)
        
        print('\nOutput:')
        print('%10s  %25s  %25s' % ('', 'rad/km/...','deg'))
        print()
        print('%10s  %25.15f' % ('Req: ', Req) )
        print('%10s  %25.15f' % ('RpolEq: ', RpolEq) )
        print()
        print('%10s  %25.15f  %25.15f' % ('u: ', u, u*_ac.r2d) )
        print('%10s  %25.15f' % ('RsinPhi: ', RsinPhi) )
        print('%10s  %25.15f' % ('RcosPhi: ', RcosPhi) )
        print()
        print('%10s  %25.15f' % ('sinHp: ', sinHp) )
        print('%10s  %25.15f  %25.15f' % ('N: ', N, N*_ac.r2d) )
        print()
        print('%10s  %25.15f  %25.15f' % ('lon_tc: ', lon_tc, lon_tc*_ac.r2d) )
        print('%10s  %25.15f  %25.15f' % ('lat_tc: ', lat_tc, lat_tc*_ac.r2d) )
        print('%10s  %25.15f  %25.15f' % ('rad_tc: ', rad_tc, rad_tc*_ac.r2d) )
        print()
        
    return lon_tc,lat_tc,rad_tc


# Test code:
if __name__ == '__main__':
    lon_tc,lat_tc,rad_tc = geoc2topoc_ecl(0.0,0.0, _ac.au,_ac.au/1000, 23*_ac.d2r,0.0, 52*_ac.d2r,0.0, debug=True)
    print(lon_tc,lat_tc,rad_tc)  # CHECK: _ac.au/1000 doesn't make sense - need apparent diameter!
    
