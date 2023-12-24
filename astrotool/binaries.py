#!/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: EUPL-1.2

#  Copyright (c) 2019-2023  Marc van der Sluys - marc.vandersluys.nl
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


"""Functions on binary stars for AstroTool."""


# Allow relative imports from __main__() when running this file (PEP 366):
if __name__ == '__main__' and __package__ is None:
    __package__ = 'astrotool'

# Modules:
import numpy as _np
import astroconst as _ac


def orb_a_from_p(m1,m2, p_orb):
    """Compute the orbital separation from the masses and orbital period using Kepler's law.
    
    Parameters:
      m1 (float):     Mass of object 1 (Mo).
      m2 (float):     Mass of object 2 (Mo).
      p_orb (float):  Orbital period (days).
    
    Returns:
      (float):  Orbital separation (Ro).
    """
    
    Mtot = (m1+m2) * _ac.m_sun
    a_orb = ((_ac.g*Mtot)/(4*_ac.pi**2))**(1/3) * (p_orb*_ac.day)**(2/3) / _ac.r_sun
    
    return a_orb


def orb_p_from_a(m1,m2, a_orb):
    """Compute the orbital period from the masses and orbital separation using Kepler's law.
    
    Parameters:
      m1 (float):     Mass of object 1 (Mo).
      m2 (float):     Mass of object 2 (Mo).
      a_orb (float):  Orbital separation (Ro).
    
    Returns:
      (float):  Orbital period (days).
    """
    
    Mtot = (m1+m2) * _ac.m_sun
    p_orb = _np.sqrt(4*_ac.pi**2/(_ac.g*Mtot)) * (a_orb*_ac.r_sun)**1.5 / _ac.day
    
    return p_orb


def orb_en_from_a(m1,m2, a_orb):
    """Compute the orbital energy from the masses and orbital separation.
    
    Parameters:
      m1 (float):     Mass of object 1 (Mo).
      m2 (float):     Mass of object 2 (Mo).
      a_orb (float):  Orbital separation (Ro).
    
    Returns:
      (float):  Orbital energy (absolute value - J).
    """
    
    return _ac.g * m1*m2*_ac.m_sun**2 / (2*a_orb*_ac.r_sun)


def orb_a_from_en(m1,m2, e_orb):
    """Compute the orbital separation from the masses and orbital energy.
    
    Parameters:
      m1 (float):     Mass of object 1 (Mo).
      m2 (float):     Mass of object 2 (Mo).
      e_orb (float):  Orbital energy (J).
    
    Returns:
      (float):  Orbital separation (Ro).
    """
    
    a_orb = (_ac.g * m1*m2*_ac.m_sun**2) / (2 * e_orb)  / _ac.r_sun  # in Ro
    
    return a_orb


def orb_am_from_a(m1,m2, a_orb):
    """Compute the orbital angular momentum from the masses and orbital separation.
    
    Parameters:
      m1 (float):     Mass of object 1 (Mo).
      m2 (float):     Mass of object 2 (Mo).
      a_orb (float):  Orbital separation (Ro).
    
    Returns:
      (float):  Orbital angular momentum (J s).
    """
    
    Mtot = (m1+m2) * _ac.m_sun
    l_orb = m1*m2*_ac.m_sun**2 * _np.sqrt(_ac.g*a_orb*_ac.r_sun/Mtot)
    
    return l_orb


def orb_a_from_am(m1,m2, l_orb):
    """Compute the orbital separation from the masses and orbital angular momentum.
    
    Parameters:
      m1 (float):     Mass of object 1 (Mo).
      m2 (float):     Mass of object 2 (Mo).
      l_orb (float):  Orbital angular momentum (J s).
    
    Returns:
      (float):  Orbital separation (Ro).
    """
    
    Mtot = (m1+m2) * _ac.m_sun
    a_orb = Mtot/_ac.g * _np.square(l_orb / (m1*m2*_ac.m_sun**2)) / _ac.r_sun
    
    return a_orb


# Test code:
if __name__ == '__main__':
    pass
    
