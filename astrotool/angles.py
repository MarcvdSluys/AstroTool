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


"""Angle functions for AstroTool."""


# Allow relative imports from __main__() when running this file (PEP 366):
if __name__ == '__main__' and __package__ is None:
    __package__ = 'astrotool'

# Modules:
import numpy as _np
import astroconst as _ac


def map_pi_pi(ang):
    """Map an (array of) angle(s) to lie between -PI and +PI.
    
    Parameters:
      ang (float):  (Array of) angle(s).
    
    Returns:
      (float):  (Array) of angles between -PI and +PI.
    """
    
    return _np.mod(ang+_ac.pi, _ac.pi2)-_ac.pi


def angle_sum(ang1, ang2):
    """Return the sum of two angles, between -pi and +pi.
    
    Parameters:
      ang1 (float):  Angle 1 (rad).
      ang2 (float):  Angle 2 (rad).
    
    Returns:
      (float):  Sum of angle 1 and angle 2 (rad).
    """
    
    asum = map_pi_pi(ang1 + ang2)
    
    return asum


def angle_diff(ang1, ang2):
    """Return the difference between two angles, between -pi and +pi.
    
    Parameters:
      ang1 (float):  Angle 1 (rad).
      ang2 (float):  Angle 2 (rad).
    
    Returns:
      (float):  Difference between angle 1 and angle 2 (rad).
    """
    
    diff = map_pi_pi(ang1 - ang2)
    
    return diff


def position_angle(lon1,lat1, lon2,lat2):
    """Calculate the position angle of object 2 with respect to object 1, counterclockwise from the north.
    
    Parameters:
      lon1 (float):  Longitude of object 1 (rad).
      lat1 (float):  Latitude of object 1 (rad).
    
      lon2 (float):  Longitude of object 2 (rad).
      lat2 (float):  Latitude of object 2 (rad).
    
    Returns:
      (float):  Position angle (rad).
    
    Note:
      You may want to swap object1 and object2 for e.g. right ascension, which increases in the other direction.
    """
    
    dlon = lon2-lon1
    pa   = _np.arctan2( _np.sin(dlon),  _np.cos(lat1)*_np.tan(lat2) - _np.sin(lat1)*_np.cos(dlon) )
    
    return pa
    

# Test code:
if __name__ == '__main__':
    pass
    
