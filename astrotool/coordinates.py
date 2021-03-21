#!/bin/env python

#  Copyright (c) 2019-2021  Marc van der Sluys - marc.vandersluys.nl
#   
#  This file is part of the AstroTool Python package,
#  see: http://astro.ru.nl/~sluys/AstroTool/
#   
#  AstroTool has been developped by Marc van der Sluys of the Department of Astrophysics at the Radboud
#  University Nijmegen, the Netherlands and the departmen of Sustainable energy of the HAN University of
#  applied sciences in Arnhem, the Netherlands.
#   
#  This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#  
#  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
#  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License along with this code.  If not, see
#  <http://www.gnu.org/licenses/>.


"""Coordinate transformations and related functions for AstroTool."""


# Modules:
import numpy.core as np
from astrotool.constants import pi2, r2d,as2r, earth_rad,AU
import astrotool.date_time as atdt



def obliquity(jd):
    """Compute the obliquity of the ecliptic in radians from the JD(E).
    
    Parameters:
      jd (double):  Julian day (days).
    
    Returns:
      double:  eps: Obliquity of the ecliptic (rad).
    
    References:
      - Seidelman 1992, Eq. 3.222-1.

    """
    
    tJC = atdt.jd2tjc(jd)  # Time in Julian centuries since J2000.0
    eps = 0.409092804 - 2.269655e-4*tJC - 2.86e-9*tJC**2 + 8.78967e-9*tJC**3  # Obliquity of the ecliptic (rad)
    
    return eps


def eq2ecl(ra,dec, eps):
    """Convert equatorial coordinates to ecliptical.

    Parameters:
      ra (double):   Right ascension (rad).
      dec (double):  Declination (rad).
      eps (double):  Obliquity of the ecliptic (rad).
    
    Returns:
      tuple (double,double): tuple containing (lon, lat):
    
        - lon (double):  Ecliptical longitude (rad).
        - lat (double):  Ecliptical latitude (rad).

    """
    
    lon = np.arctan2( np.sin(ra)  * np.cos(eps) + np.tan(dec) * np.sin(eps),  np.cos(ra) ) % pi2
    lat =  np.arcsin( np.sin(dec) * np.cos(eps) - np.cos(dec) * np.sin(eps) * np.sin(ra) )
    
    return lon,lat


def ecl2eq(lon,lat, eps):
    """Convert (geocentric) spherical ecliptical coordinates to spherical equatorial coordinates.
    
    Parameters:
      lon (double):  Ecliptical longitude (rad).
      lat (double):  Ecliptical latitude (rad).
      eps (double):  Obliquity of the ecliptic (rad).
    
    Returns:
      tuple (double,double): tuple containing (ra, dec):
    
        - ra (double):   Right ascension (rad).
        - dec (double):  Declination (rad).
    
    References:
      - `Explanatory Supplement to the Astronomical Almanac 3rd Ed,
        Eq.14.43 <https://aa.usno.navy.mil/publications/docs/exp_supp.php>`_

    """
    
    ra  = np.arctan2( np.sin(lon) * np.cos(eps)  -  np.tan(lat) * np.sin(eps),  np.cos(lon) ) % pi2
    dec =  np.arcsin( np.sin(lat) * np.cos(eps)  +  np.cos(lat) * np.sin(eps) * np.sin(lon) )
    
    return ra,dec


def par2horiz(ha,dec, phi):
    """Convert parallactic coordinates to horizontal.
    
    Parameters:
      ha (double):   Hour angle (rad).
      dec (double):  Declination (rad).
      phi (double):  Geographical latitude (rad, N>0).
    
    Returns:
      tuple (double,double): tuple containing (az, alt):
    
        - az (double):   Azimuth (rad, S=0).
        - alt (double):  Altitude (rad).
    
    """
    
    az  = np.arctan2( np.sin(ha),   np.cos(ha) * np.sin(phi) - np.tan(dec) * np.cos(phi) ) % pi2
    alt = np.arcsin(  np.sin(dec) * np.sin(phi) + np.cos(ha) * np.cos(dec) * np.cos(phi) )
    
    return az,alt


def proper_motion(jd_start,jd_target, ra,dec, pma,pmd):
    """Compute the proper motion from jd_start to jd_target for the given positions and proper motions.
    
    Parameters:
      jd_start (double):   Julian day of the initial epoch (days).
      jd_target (double):  Julian day of the target epoch (days).
    
      ra (double):         Right ascension (numpy array, rad).
      dec (double):        Declination (numpy array, rad).
    
      pma (double):        Proper motion in right ascension (numpy array, rad/yr).
      pmd (double):        Proper motion in declination (numpy array, rad/yr).
    
    Returns:
      tuple (double,double): tuple containing (ra_target, dec_target):
    
        - ra_target (double):   Right ascension for the target epoch (rad).
        - dec_target (double):  Declination for the target epoch (rad).
    
    """
    
    dtYr   = (jd_target - jd_start)/365.25
    ra_target  = ra  + pma*dtYr / np.cos(dec)
    dec_target = dec + pmd*dtYr
    
    return ra_target,dec_target


def precess_from_2000(jd, ra,dec):
    """Compute precession in equatorial coordinates from J2000 to that of the specified JD.
    
    J2000 is the equinox of many catalogues, including the Hipparcos one.
    
    Parameters:
      jd (double):   Julian day (days).
      ra (double):   Right ascension (rad).
      dec (double):  Declination (rad).
    
    Returns:
      tuple (double,double): tuple containing (ra_new, dec_new):
    
        - ra_new (double):   Right ascension for the target equinox (rad).
        - dec_new (double):  Declination for the target equinox (rad).
    
    """
    
    tJC  = atdt.jd2tjc(jd)  # Time in Julian centuries since J2000.0
    tJC2 = tJC**2
    tJC3 = tJC*tJC2
    
    zeta  = (2306.2181*tJC + 0.30188*tJC2 + 0.017998*tJC3)*as2r  # Convert the result from arcseconds to radians
    z     = (2306.2181*tJC + 1.09468*tJC2 + 0.018203*tJC3)*as2r
    theta = (2004.3109*tJC - 0.42665*tJC2 - 0.041833*tJC3)*as2r
    
    ra_new  = (np.arctan2( np.sin(ra + zeta) * np.cos(dec),  np.cos(ra + zeta) * np.cos(theta) * np.cos(dec) - np.sin(theta) * np.sin(dec) ) + z) % pi2
    dec_new = np.arcsin( np.cos(ra + zeta) * np.sin(theta) * np.cos(dec)  +  np.cos(theta) * np.sin(dec) )
    
    return ra_new,dec_new


def geoc2topoc_ecl(lon_gc,lat_gc, dist_gc,rad_gc, eps,lst, lat_obs,ele_obs=0, debug=False):
    """Convert spherical ecliptical coordinates from the geocentric to the topocentric system.
    
    Parameters:
      lon_gc (double):   Geocentric ecliptic longitude (rad).
      lat_gc (double):   Geocentric ecliptic latitude (rad).
      dist_gc (double):  Geocentric distance (AU).
      rad_gc (double):   Geocentric semi-diameter (rad).
      
      eps (double):      Obliquity of the ecliptic (rad).
      lst (double):      Local sidereal time (rad).
      
      lat_obs (double):  Geographical latitude of the observer (rad).
      ele_obs (double):  Altitude/elevation of the observer above sea level (metres, optional, default value = 0).
      
      debug (bool):      Print debug output (True/False, optional, default value = True).
    
    Returns:
      tuple (double,double,double): tuple containing (lon_tc, lat_tc, rad_tc):
    
        - lon_tc (double):  Topocentric ecliptic longitude (rad).
        - lat_tc (double):  Topocentric ecliptic latitude (rad).
        - rad_tc (double):  Topocentric semi-diameter (rad).
    
    """
    
    # Meeus, Ch.11, p.82:
    Req = earth_rad*1000      # Equatorial radius of the Earth in metres (same units as the elevation)
    #                         (http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf)
    RpolEq = 0.996647189335   # Rpol/Req = 1-f: flattening of the Earth - WGS84 ellipsoid 
    
    u  = np.arctan(RpolEq*np.tan(lat_obs))
    RsinPhi = RpolEq*np.sin(u) + ele_obs/Req * np.sin(lat_obs)
    RcosPhi = np.cos(u)        + ele_obs/Req * np.cos(lat_obs)
    
    sinHp = np.sin(earth_rad/AU)/(dist_gc/AU)  # Sine of the horizontal parallax, Meeus, Eq. 40.1
    
    # Meeus, Ch.40, p.282:
    N  = np.cos(lon_gc)*np.cos(lat_gc) - RcosPhi*sinHp*np.cos(lst)
    
    lon_tc = np.arctan2( np.sin(lon_gc)*np.cos(lat_gc) - sinHp*(RsinPhi*np.sin(eps) + RcosPhi*np.cos(eps)*np.sin(lst)), N ) % pi2  # Topocentric longitude
    lat_tc = np.arctan((np.cos(lon_tc)*(np.sin(lat_gc) - sinHp*(RsinPhi*np.cos(eps) - RcosPhi*np.sin(eps)*np.sin(lst))))/N)        # Topocentric latitude
    rad_tc = np.arcsin(np.cos(lon_tc)*np.cos(lat_tc)*np.sin(rad_gc)/N)                                                             # Topocentric semi-diameter
    
    # print(dist_gc, dist_gc*rad_gc/rad_tc)
    
    if(debug):
        print()
        print('geoc2topoc_ecl():')
        print('%10s  %25s  %25s' % ('', 'rad/km/...','deg'))
        print()
        print('%10s  %25.15f' % ('Req: ', Req) )
        print('%10s  %25.15f' % ('RpolEq: ', RpolEq) )
        print()
        print('%10s  %25.15f  %25.15f' % ('u: ', u, u*r2d) )
        print('%10s  %25.15f' % ('RsinPhi: ', RsinPhi) )
        print('%10s  %25.15f' % ('RcosPhi: ', RcosPhi) )
        print()
        print('%10s  %25.15f' % ('sinHp: ', sinHp) )
        print('%10s  %25.15f  %25.15f' % ('N: ', N, N*r2d) )
        print()
        print('%10s  %25.15f  %25.15f' % ('lon_tc: ', lon_tc, lon_tc*r2d) )
        print('%10s  %25.15f  %25.15f' % ('lat_tc: ', lat_tc, lat_tc*r2d) )
        print('%10s  %25.15f  %25.15f' % ('rad_tc: ', rad_tc, rad_tc*r2d) )
        print()
        
    return lon_tc,lat_tc,rad_tc


# Test code:
if(__name__ == "__main__"):
    from constants import d2r
    print(geoc2topoc_ecl(0.0,0.0, AU,AU, 23*d2r,0.0, 52*d2r,0.0, debug=True))
    
