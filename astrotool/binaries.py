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


"""Functions on binary stars for AstroTool."""


# Allow relative imports from __main__() when running this file (PEP 366):
if __name__ == '__main__' and __package__ is None:
    __package__ = 'astrotool'

# Modules:
import numpy as _np
import astroconst as _ac


def orb_a_from_p(m1,m2, Porb):
    """Compute the orbital separation from the masses and orbital period using Kepler's law.
    
    Parameters:
      m1 (float):    Mass of star 1 (Mo).
      m2 (float):    Mass of star 2 (Mo).
      Porb (float):  Orbital period (days).
    
    Returns:
      (float):  Orbital separation (Ro).
    """
    
    Mtot = (m1+m2) * _ac.m_sun
    a_orb = ((_ac.g*Mtot)/(4*_ac.pi**2))**(1/3) * (Porb*_ac.day)**(2/3) / _ac.r_sun
    
    return a_orb


def orb_p_from_a(m1,m2, a_orb):
    """Compute the orbital period from the masses and orbital separation using Kepler's law.
    
    Parameters:
      m1 (float):     Mass of star 1 (Mo).
      m2 (float):     Mass of star 2 (Mo).
      a_orb (float):  Orbital separation (Ro).
    
    Returns:
      (float):  Orbital period (days).
    """
    
    Mtot = (m1+m2) * _ac.m_sun
    Porb = _np.sqrt(4*_ac.pi**2/(_ac.g*Mtot)) * (a_orb*_ac.r_sun)**1.5 / _ac.day
    
    return Porb


def orb_en_from_a(m1,m2, a_orb):
    """Compute the orbital energy from the masses and orbital separation.
    
    Parameters:
      m1 (float):     Mass of star 1 (Mo).
      m2 (float):     Mass of star 2 (Mo).
      a_orb (float):  Orbital separation (Ro).
    
    Returns:
      (float):  Orbital energy (absolute value - J).
    """
    
    return _ac.g * m1*m2*_ac.m_sun**2 / (2*a_orb*_ac.r_sun)


def orb_a_from_en(m1,m2, e_orb):
    """Compute the orbital separation from the masses and orbital energy.
    
    Parameters:
      m1 (float):     Mass of star 1 (Mo).
      m2 (float):     Mass of star 2 (Mo).
      e_orb (float):  Orbital energy (J).
    
    Returns:
      (float):  Orbital separation (Ro).
    """
    
    a_orb = (_ac.g * m1*m2*_ac.m_sun**2) / (2 * e_orb)  / _ac.r_sun  # in Ro
    
    return a_orb


def orb_am_from_a(m1,m2, a_orb):
    """Compute the orbital angular momentum from the masses and orbital separation.
    
    Parameters:
      m1 (float):     Mass of star 1 (Mo).
      m2 (float):     Mass of star 2 (Mo).
      a_orb (float):  Orbital separation (Ro).
    
    Returns:
      (float):  Orbital angular momentum (SI: J s).
    """
    
    Mtot = (m1+m2) * _ac.m_sun
    Jorb = m1*m2*_ac.m_sun**2 * _np.sqrt(_ac.g*a_orb*_ac.r_sun/Mtot)
    
    return Jorb


def orb_a_from_am(m1,m2, Jorb):
    """Compute the orbital separation from the masses and orbital angular momentum.
    
    Parameters:
      m1 (float):    Mass of star 1 (Mo).
      m2 (float):    Mass of star 2 (Mo).
      Jorb (float):  Orbital angular momentum (SI: J s).
    
    Returns:
      (float):  Orbital separation (Ro).
    """
    
    Mtot = (m1+m2) * _ac.m_sun
    a_orb = Mtot/_ac.g * _np.square(Jorb / (m1*m2*_ac.m_sun**2)) / _ac.r_sun
    
    return a_orb


def roche_lobe_from_a_pac(m1,m2, a_orb):
    """Compute the Roche lobe for star 1 from the masses and orbital separation using Paczynski (1967).
    
    Parameters:
      m1 (float):     Mass of star 1, for which the Roche lobe should be computed (arbitrary).
      m2 (float):     Mass of star 2 (same as m1).
      a_orb (float):  Orbital separation (Ro).
    
    Returns:
      (float):  Roche lobe of the star with m1 (same as a_orb).
    
    Note: fixed typo from the original paper.
    """
    
    return 2/(3**(4/3)) * a_orb * (m1/(m1+m2))**(1/3)


def roche_lobe_from_a_egg(m1,m2, a_orb):
    """Compute the Roche lobe for star 1 from the masses and orbital separation using Eggleton (1983).
    
    Parameters:
      m1 (float):     Mass of star 1, for which the Roche lobe should be computed (arbitrary).
      m2 (float):     Mass of star 2 (same as m1).
      a_orb (float):  Orbital separation (Ro).
    
    Returns:
      (float):  Roche lobe of the star with m1 (same as a_orb).
    """
    
    q1 = m1/m2
    q123 = q1**(2/3)
    r_rl = a_orb * 0.49*q123 / (0.6*q123 + _np.log(1+q1**(1/3)))  # Note: ln()
    
    return r_rl


def roche_lobe_from_a_egg_simpl(m1,m2, a_orb):
    """Compute the Roche lobe for star 1 from the masses and orbital separation using simplified Eggleton (2006).
    
    Parameters:
      m1 (float):     Mass of star 1, for which the Roche lobe should be computed (arbitrary).
      m2 (float):     Mass of star 2 (same as m1).
      a_orb (float):  Orbital separation (Ro).
    
    Returns:
      (float):  Roche lobe of the star with m1 (same as a_orb).
    """
    
    q1 = m1/m2
    
    return a_orb * 0.44*q1**0.33 / (1 + q1)**0.2


def roche_lobe_accel_l1(r2, m1,m2, aorb):
    """Return the acceleration (the gradient of the potential = F/m) for a distance r2 "in front of" star 2
    (as seen from the centre of mass) on the line connecting the stars with masses m1 and m2 of the binary
    with orbital separation aorb.  This is useful to find the position of the first Lagrangian point L1.
    
    Parameters:
      r2 (float):    Distance "before" star 2, as seen from the centre of mass, in terms of aorb.
      m1 (float):    Mass of star 1.
      m2 (float):    Mass of star 2.
      aorb (float):  Orbital separation of the two stars.
    
    Returns:
      (float):  acceleration at the given point.
    """
    
    mtot = m1+m2
    a2 = m1/mtot * aorb
    accel = m1/(aorb-r2)**2 - m2/r2**2 - (a2-r2)*mtot/aorb**3
    return accel


def roche_lobe_accel_l2(r2, m1,m2, aorb):
    """Return the acceleration (the gradient of the potential = F/m) for a distance r2 "behind" star 2 (as
    seen from the centre of mass) on the line connecting the stars with masses m1 and m2 of the binary with
    orbital separation aorb.  This is useful to find the position of the second Lagrangian point L2.
    
    Parameters:
      r2 (float):    Distance "behind" star 2, as seen from the centre of mass, in terms of aorb.
      m1 (float):    Mass of star 1.
      m2 (float):    Mass of star 2.
      aorb (float):  Orbital separation of the two stars.
    
    Returns:
      (float):  acceleration at the given point.
    """
    
    mtot = m1+m2
    a2 = m1/mtot * aorb
    accel = -m1/(aorb+r2)**2 - m2/r2**2 + (a2+r2)*mtot/aorb**3
    return accel


def roche_lobe_accel_l3(r1, m1,m2, aorb):
    """Return the acceleration (the gradient of the potential = F/m) for a distance r3 "behind" star 1 (as
    seen from the centre of mass) on the line connecting the stars with masses m1 and m2 of the binary with
    orbital separation aorb.  This is useful to find the position of the second Lagrangian point L3.
    
    Parameters:
      r1 (float):    Distance "behind" star 1, as seen from the centre of mass, in terms of aorb.
      m1 (float):    Mass of star 1.
      m2 (float):    Mass of star 2.
      aorb (float):  Orbital separation of the two stars.
    
    Returns:
      (float):  acceleration at the given point.
    """
    
    mtot = m1+m2
    a1 = m2/mtot * aorb
    accel = -m1/(aorb-r1)**2 - m2/(2*aorb-r1)**2 + (a1+aorb-r1)*mtot/aorb**3
    return accel


def l2_from_q2(q2):
    """Return an estimate of the distance between the second Lagrangian point and the centre of mass expressed
    in a_orb, from the mass ratio q2=M2/M1.
    
    Parameters:
      q2 (float):  The mass ratio 0 < M2/M1 <= 1.
    
    Returns:
      (float):  Distance between centre of mass and L2 in units of the binary orbital separation (a_orb).
    
    Notes:
      - This model is a fit to more detailed calculations, with rounded off coefficients.  The fit is
        cut into two pieces, below and above q2=0.19:
        - q2<0.19:  function is an "ellipse" + linear term;
        - q2>=0.19: function is a third-order polynomial;
        - the function is quite continuous at q2=0.19, though its derivative is not.
      - Accuracy:
        - Mean/med. |relative difference|:  dl2/l2 = 0.059% / 0.035     (orig fit: 0.025% / 0.024%);
        - Max. |relative difference|:       dl2/l2 = 0.306% at q2=1e-6  (orig fit: 0.144% at q2=1).
    """
    
    q2 = _np.asarray(_np.copy(q2), dtype=float)  # Copy and typecast to numpy.ndarray - ensure float!
    
    scalar_input = False
    if q2.ndim == 0:
        q2 = q2[None]  # Makes q2 a 1D array.
        scalar_input = True
    
    l2 = _np.zeros_like(q2, dtype=float)  # Ensure float!
    # l2[q2<0.19]  = 0.997205 + 0.121484  * q2[q2<0.19]  + 0.251220 * _np.power(1.00813 - _np.power((0.202059-q2[q2<0.19])/0.201651, 3.99641), 0.315883)  # "Ellipse + linear" - orig fit
    # l2[q2>=0.19] = 1.26849  + 0.0717643 * q2[q2>=0.19] - 0.296485 * _np.square(q2[q2>=0.19]) + 0.156368 * _np.power(q2[q2>=0.19], 3)  # 3rd-order polynomial - orig fit
    l2[q2<0.19]  = 0.9972 + 0.1215  * q2[q2<0.19]  + 0.2512 * _np.power(1.0081 - _np.power((0.202059-q2[q2<0.19])/0.201651, 3.996), 0.316)  # "Ellipse + linear" - step 5: ~2x worse
    l2[q2>=0.19] = 1.268 + 0.072 * q2[q2>=0.19] - 0.2965 * _np.square(q2[q2>=0.19]) + 0.1564 * _np.power(q2[q2>=0.19], 3)  # 3rd-order polynomial - step 6: slightly better (after step 5)
    
    if scalar_input:
        return float(_np.squeeze(l2))  # Array -> scalar (float)
    
    return l2


def gw_am_loss_from_a(m1,m2, aorb):
    """Return the angular-momentum loss due to GWs for a binary with given masses and orbital separation.
    
    Parameters:
      m1 (float):    Mass of star 1 (Mo).
      m2 (float):    Mass of star 2 (Mo).
      aorb (float):  Orbital separation (Ro).
    
    Returns:
      (float):  Angular-momentum loss dJ/dt (SI: J).
    """
    
    mtot = m1 + m2
    dJdt = -32/5 * _ac.G**(7/2) / _ac.c**5 * \
        _np.square(m1 * m2 * _ac.Mo**2) * _np.sqrt(mtot*_ac.Mo) / _np.power(aorb*_ac.Ro, 7/2)
    
    return dJdt
    

def gw_am_loss_from_P(m1,m2, Porb):
    """Return the angular-momentum loss due to GWs for a binary with given masses and orbital separation.
    
    Parameters:
      m1 (float):    Mass of star 1 (Mo).
      m2 (float):    Mass of star 2 (Mo).
      Porb (float):  Orbital period (days).
    
    Returns:
      (float):  Angular-momentum loss dJ/dt (SI: J).
    """
    
    mtot = m1 + m2
    dJdt = -32/5 * _ac.G**(7/3) / _ac.c**5 * \
        _np.square(m1 * m2 * _ac.Mo**2 / _np.power(mtot*_ac.Mo, 1/3)) * _np.power(_ac.pi2 / (Porb*_ac.day), 7/3)
    
    return dJdt
    

def gw_merger_time_from_a(m1,m2, aorb):
    """Return the merger time due to GWs for a binary with given masses and orbital separation.
    
    Parameters:
      m1 (float):    Mass of star 1 (Mo).
      m2 (float):    Mass of star 2 (Mo).
      aorb (float):  Orbital separation (Ro).
    
    Returns:
      (float):  Merger time (yr).
    """
    
    t_merge = 5/256 * _ac.c**5 / _ac.G**3 * _np.power(aorb * _ac.Rsun, 4) / (m1 * m2 * (m1+m2) * _ac.Msun**3)
    
    return t_merge / _ac.year
    

def gw_merger_time_from_P(m1,m2, Porb):
    """Return the merger time due to GWs for a binary with given masses and orbital period.
    
    Parameters:
      m1 (float):    Mass of star 1 (Mo).
      m2 (float):    Mass of star 2 (Mo).
      Porb (float):  Orbital period (days).
    
    Returns:
      (float):  Merger time (yr).
    """
    
    aorb = orb_a_from_p(m1,m2, Porb)
    t_merge = 5/256 * _ac.c**5 / _ac.G**3 * _np.power(aorb * _ac.Rsun, 4) / (m1 * m2 * (m1+m2) * _ac.Msun**3)
    
    return t_merge / _ac.year
    

# Test code:
if __name__ == '__main__':
    pass
    
