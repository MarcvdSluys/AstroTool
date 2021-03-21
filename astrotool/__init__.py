#!/bin/env python

#  Copyright (c) 2020-2021  Marc van der Sluys - marc.vandersluys.nl
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


"""
AstroTool package
=================

AstroTool is a Python package to do basic astronomical calculations in Python or on the command line.  The
package can be used under the conditions of the GPLv3 licence.  These pages contain the API documentation.
For more information on the Python package, licence, source code and data files, see the `AstroTool
homepage <http://astro.ru.nl/~sluys/AstroTool/>`_.

"""


name = "astrotool"

from .constants import *
from .coordinates import *
from .date_time import *

