#!/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: EUPL-1.2

#  Copyright (c) 2019-2022  Marc van der Sluys - marc.vandersluys.nl
#   
#  This file is part of the AstroTool Python package,
#  see: http://astro.ru.nl/~sluys/AstroTool/
#   
#  AstroTool has been developed by Marc van der Sluys of the Department of Astrophysics at the Radboud
#  University Nijmegen, the Netherlands and the department of Sustainable energy of the HAN University of
#  applied sciences in Arnhem, the Netherlands.
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


"""Definition of constants for AstroTool."""

# Modules:
import numpy as np

pi   = np.pi;           """pi"""
pi2  = 2*pi;            """2 pi"""
pio2 = pi/2;            """pi/2"""

r2d  = np.degrees(1);   """Radians to degrees"""
d2r  = np.radians(1);   """Degrees to radians"""

h2r   = d2r*15;         """Hours to radians"""
as2r  = d2r/3.6e3;      """Arcseconds to radians"""
mas2r = as2r/1000.0;    """Milliarcseconds to radians"""

jd1820 = 2385801;       """JD in 1820  (for Delta-T fit)"""
jd1900 = 2415021;       """JD in 1900"""
jd2000 = 2451545;       """JD in 2000.0"""

earth_rad = 6378.1366;  """Earth radius in km"""
moon_rad  = 1737.5;     """Moon radius in km"""

AU = 1.495978707e8;     """Astronomical unit in km"""


# Test code:
if __name__ == '__main__':
    print(pi, r2d, jd2000, earth_rad, AU)
    
