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


"""Functions on simple cosmology for AstroTool."""


# Allow relative imports from __main__() when running this file (PEP 366):
if __name__ == '__main__' and __package__ is None:
    __package__ = 'astrotool'

# Modules:
import numpy as _np
import astroconst as _ac
import scipy.integrate as _sp_int


def Hubble_parameter_from_z(z, Omega_m=0.3, Omega_k=0, Omega_L=0.7):
    """Compute the dimensionless Hubble parameter.
    
    Parameters:
      z (float):        Cosmological redshift.
      Omega_m (float):  Mass-density parameter ρ_m/ρ_m,crit.  Optional; defaults to 0.3.
      Omega_k (float):  Curvature-density parameter.  Optional; defaults to 0: flat universe.
      Omega_L (float):  Density parameter for the cosmological constant.  Optional; defaults to 0.7.
    
    Returns:
      (float):  Dimensionless Hubble parameter.
    """
    
    HP = _np.sqrt(Omega_m * _np.power(1 + z, 3) + Omega_k * _np.square(1 + z) + Omega_L)
    
    return HP


def _comoving_distance_from_z_scalar(z, H0=70, Omega_m=0.3, Omega_k=0, Omega_L=0.7):
    """Compute the comoving radial distance in Mpc from the cosmological redshift.
    
    Parameters:
      z (float):        Cosmological redshift.
      H0 (float):       Hubble constant in km/s/Mpc.
      Omega_m (float):  Mass-density parameter ρ_m/ρ_m,crit.  Optional; defaults to 0.3.
      Omega_k (float):  Curvature-density parameter.  Optional; defaults to 0: flat universe.
      Omega_L (float):  Density parameter for the cosmological constant.  Optional; defaults to 0.7.
    
    Returns:
      (float):  comoving radial distance (Mpc).
    """
    
    integral, _ = _sp_int.quad(lambda zp: 1.0 / Hubble_parameter_from_z(zp, Omega_m, Omega_k, Omega_L), 0, z)  # d_C(z) = d_H ∫_0^z E(z') dz'
    integral *= (_ac.c/1000 / H0)
    
    return integral  


def _angular_diameter_distance_from_z_scalar(z):
    """Compute the angular-diameter distance in Mpc from the cosmological redshift.
    
    Parameters:
      z (float):  Cosmological redshift.
    
    Returns:
      (float):  angular-diameter distance (Mpc).
    """
    
    dist = _comoving_distance_from_z_scalar(z) / (1 + z)
    
    return dist


def _luminosity_distance_from_z_scalar(z):
    """Compute the luminosity distance in Mpc from the cosmological redshift.
    
    Parameters:
      z (float):  Cosmological redshift.
    
    Returns:
      (float):  angular-diameter distance (Mpc).
    """
    
    dist = _comoving_distance_from_z_scalar(z) * (1 + z)
    
    return dist


def _lookback_time_from_z_scalar(z, H0=70, Omega_m=0.3, Omega_k=0, Omega_L=0.7):
    """Compute the lookback time in Gyr from the cosmological redshift.
    
    Parameters:
      z (float):        Cosmological redshift.
      H0 (float):       Hubble constant in km/s/Mpc.
      Omega_m (float):  Mass-density parameter ρ_m/ρ_m,crit.  Optional; defaults to 0.3.
      Omega_k (float):  Curvature-density parameter.  Optional; defaults to 0: flat universe.
      Omega_L (float):  Density parameter for the cosmological constant.  Optional; defaults to 0.7.
    
    Returns:
      (float):  lookback time (Gyr).
    """
    
    H0_SI = H0 * 1e3 / _ac.Mpc   # H0 in km/s/Mpc -> 1/s
    t_H = 1 / H0_SI / _ac.Gyr    # Hubble time in s -> Gyr
    integral, _ = _sp_int.quad(lambda zp: 1.0 / ((1 + zp) * Hubble_parameter_from_z(zp, Omega_m, Omega_k, Omega_L)), 0, z)
    time = t_H * integral
    
    return time


def _lookback_distance_from_z_scalar(z, H0=70, Omega_m=0.3, Omega_k=0, Omega_L=0.7):
    """Compute the lookback distance (light-travel distance) in Mpc from the cosmological redshift.
    
    Parameters:
      z (float):        Cosmological redshift.
      H0 (float):       Hubble constant in km/s/Mpc.
      Omega_m (float):  Mass-density parameter ρ_m/ρ_m,crit.  Optional; defaults to 0.3.
      Omega_k (float):  Curvature-density parameter.  Optional; defaults to 0: flat universe.
      Omega_L (float):  Density parameter for the cosmological constant.  Optional; defaults to 0.7.
    
    Returns:
      (float):  lookback distance (Mpc).
    
    Note:
      This is the lookback time converted to Mpc.
    """
    
    t_lb = _lookback_time_from_z_scalar(z, H0, Omega_m, Omega_k, Omega_L) * _ac.Gyr  # Gyr -> s
    dist = _ac.c * t_lb / _ac.Mpc  # c t in m -> Mpc
    
    return dist


# Vectorize functions above for arrays:
comoving_distance_from_z          = _np.vectorize(_comoving_distance_from_z_scalar)
angular_diameter_distance_from_z  = _np.vectorize(_angular_diameter_distance_from_z_scalar)
luminosity_distance_from_z        = _np.vectorize(_luminosity_distance_from_z_scalar)
lookback_time_from_z              = _np.vectorize(_lookback_time_from_z_scalar)
lookback_distance_from_z          = _np.vectorize(_lookback_distance_from_z_scalar)
    


# Test code:
if __name__ == '__main__':
    pass
    
