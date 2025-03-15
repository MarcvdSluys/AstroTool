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


import argparse
from . import stars
from . import binaries

"""CLI scripts/entrypoints for AstroTool."""


def orb_a_from_p():
    """Cli wrapper for orb_a_from_p()."""
    
    # Parse command-line arguments:
    parser = argparse.ArgumentParser(description='Compute binary orbital separation from orbital period.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)  # Use capital, period, add default values
    
    # Required arguments:
    parser.add_argument('M1',   type=float, help='mass of star 1 (Mo)')
    parser.add_argument('M2',   type=float, help='mass of star 2 (Mo)')
    parser.add_argument('Porb', type=float, help='binary orbital period (days)')
    
    # Optional arguments:
    # Mutually exclusive arguments:
    verb  = parser.add_mutually_exclusive_group()
    verb.add_argument('-v', '--verbosity', action='count', default=1, help='increase output verbosity')  # Counts number of occurrences
    verb.add_argument('-q', '--quiet',     action='store_true', help='print result only')  # False by default
    
    args = parser.parse_args()
    if args.quiet: args.verbosity=0
    
    aorb = binaries.orb_a_from_p(args.M1, args.M2, args.Porb)
    Rrl1 = binaries.roche_lobe_from_a_egg(args.M1, args.M2, aorb)
    Rrl2 = binaries.roche_lobe_from_a_egg(args.M2, args.M1, aorb)
    
    if args.verbosity>0:
        print('M1    = %0.3f Mo' % (args.M1))
        print('M2    = %0.3f Mo' % (args.M2))
        print('Porb  = %0.3f Mo' % (args.Porb))
        print()
        print('q1    = %0.3f Mo' % (args.M1/args.M2))
        print('q2    = %0.3f Mo' % (args.M2/args.M1))
        print()
        print('a_orb = %0.3f Ro' % (aorb))
        print('R_Rl1 = %0.3f Ro' % (Rrl1))
        print('R_Rl2 = %0.3f Ro' % (Rrl2))
    else:
        print(aorb)
    
    return


def wd_radius_from_mass():
    """Cli wrapper for wd_radius_from_mass()."""
    
    # Parse com`mand-line arguments:
    parser = argparse.ArgumentParser(description='Compute a white-dwarf radius from its mass.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)  # Use capital, period, add default values
    
    # Required arguments:
    parser.add_argument('mass',   type=float, help='mass of the white dwarf (Mo)')
    
    parser.add_argument('-v', '--verbosity', action='count', default=0, help='increase output verbosity')  # Counts number of occurrences
    
    args = parser.parse_args()
    
    rad = stars.wd_radius_from_mass(args.mass)
    if args.verbosity>0:
        print('White-dwarf mass:   ', args.mass, 'Mo')
        print('White-dwarf radius: ', rad,       'Ro  =  ',rad/0.00915375, 'Rearth')
    else:
        print(rad)
    
    return


# Test code:
if __name__ == '__main__':
    pass
    
