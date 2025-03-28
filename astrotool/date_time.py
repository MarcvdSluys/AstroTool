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


"""Date and time functions for AstroTool."""


# Allow relative imports from __main__() when running this file (PEP 366):
if __name__ == '__main__' and __package__ is None:
    __package__ = 'astrotool'

# Modules:
import numpy as _np
import pandas as _pd
import astroconst as _ac


def julian_day(year,month,day, jd_start_greg=2299160.5):
    """Obsolescent function.  Use jd_from_date() instead."""
    _warn_obsolescent('julian_day', 'jd_from_date', rename=True)
    return jd_from_date(year,month,day, jd_start_greg)


def jd_from_date(year,month,day, jd_start_greg=2299160.5):
    """Compute the Julian Day from given year, month and day.
    
    Args:
      year (int):             Year CE.  Note that year=0 = 1 BCE, year=-1 = 2 BCE, etc.
      month (int):            Month number of year (1-12).
      day (float):            Day of month with fraction (1.0-31.999...).
    
      jd_start_greg (float):  JD of start of Gregorian calendar (optional; default=2299160.5 = 1582-10-15.0).
    
    Returns:
      float:  jd: Julian day (days).
      
    Note:
      - The JD will be in the same timezone as the date and time (use UTC for the offical JD).
      - Decimals can be used in the day to take into account the time of day other than midnight, e.g. 1.5 for
        noon on the first day of the month.
    """
    
    # Copy and typecast input to numpy.ndarrays:
    year  = _np.asarray(_np.copy(year))
    month = _np.asarray(_np.copy(month))
    day   = _np.asarray(_np.copy(day))
    
    # Jan/Feb are month 13/14 of the previous year:
    year[month  <= 2] -= 1
    month[month <= 2] += 12
    
    # JD for Julian calendar (ensure it always is an array):
    jd = _np.asarray(_np.floor(365.25*(year+4716)) + _np.floor(30.6001*(month+1)) + day - 1524.5)
    
    # Apply correction for Gregorian calendar:
    sel      = jd >= jd_start_greg                 # Select cases for Greg.cal.
    cent_1   = _np.floor(year[sel]/100.0)          # Count: (century - 1)
    jd[sel] += 2 - cent_1 + _np.floor(cent_1/4.0)  # Offset Julian-Gregorian
    
    return _np.squeeze(jd)  # Turn array back into scalar if needed


def date_time2jd(year,month,day, hour,minute,second):
    """Obsolescent function.  Use jd_from_date_time() instead."""
    _warn_obsolescent('date_time2jd', 'jd_from_date_time', rename=True, extra=True)
    return jd_from_date_time(year,month,day, hour,minute,second)


def jd_from_date_time(year,month,day, hour,minute,second, jd_start_greg=2299160.5):
    """Compute the Julian Day from given year, month, day, hour, minute and second.
    
    Args:
      year (int):             Year CE.  Note that year=0 = 1 BCE, year=-1 = 2 BCE, etc.
      month (int):            Month number of year (1-12).
      day (int):              Day of month (1-31).
    
      hour (int):             Hour of day (0-23).
      minute (int):           Minute of hour (0-59).
      second (float):         Second of minute (0.0-59.999...).
    
      jd_start_greg (float):  JD of start of Gregorian calendar (optional; default=2299160.5 = 1582-10-15.0).
    
    Returns:
      float:  jd: Julian day (days).
      
    Note:
      - The JD will be in the same timezone as the date and time (use UTC for the offical JD).
      - Decimals can be used in the second.
    """
    
    day_f = day + hour/24 + minute/1440 + second/86400    # Day with time as fraction
    return jd_from_date(year,month,day_f, jd_start_greg)


def jd_from_datetime(dtm, jd_start_greg=2299160.5):
    """Compute the Julian Day from given year, month, day, hour, minute and second.
    
    Args:
      dtm (dt.datetime):      Datetime object.
      jd_start_greg (float):  JD of start of Gregorian calendar (optional; default=2299160.5 = 1582-10-15.0).
    
    Returns:
      float:  jd: Julian day (days).
      
    Note:
      - The JD will be in the same timezone as the date and time (use UTC for the offical JD).
    """
    
    dtm = _np.asarray(_np.copy(dtm))    # Copy and typecast to numpy.ndarray
    if dtm.ndim == 0:  dtm = dtm[None]  # Makes dtm 1D.  Comment: use np.newaxis instead?
    
    dtm = _pd.to_datetime(dtm)
    day_f = dtm.day + dtm.hour/24 + dtm.minute/1440 + dtm.second/86400    # Day with time as fraction
    jd = jd_from_date(dtm.year, dtm.month, day_f, jd_start_greg)
    
    return _np.squeeze(jd)  # Arrays -> "scalars"


def jd2ymd(jd, jd_start_greg=2299160.5):
    """Obsolescent function.  Use date_from_jd() instead."""
    _warn_obsolescent('jd2ymd', 'date_from_jd', rename=True)
    return date_from_jd(jd, jd_start_greg)


def date_from_jd(jd, jd_start_greg=2299160.5):
    """Compute the calendar date from a given Julian Day.
    
    Args:
      jd (float):             Julian day (days).
      jd_start_greg (float):  JD of start of Gregorian calendar (optional; default=2299160.5 = 1582-10-15.0).
    
    Returns:
      tuple (int,int,float):  Tuple containing (year, month, day):
    
        - year (int):   Year CE.  Note that year=0 indicates 1 BCE.
        - month (int):  Month number of year (1-12).
        - day (float):  Day of month with fraction (1.0-31.999...).
    
    Note:
      - Date and time will be in the same timezone as the JD (UT for the offical JD).
      - Decimals can be returned in the day to indicate the time of day, e.g. 1.0 for midnight and 1.5 for
        noon on the first day of the month.
    """
    
    jd  = _np.asarray(_np.copy(jd))    # Copy and typecast to numpy.ndarray
    if jd.ndim == 0:  jd = jd[None]    # Create an array from a scalar if needed
    
    jd05  = _np.floor(jd+0.5)  # Julian day + 0.5 -> .0 at midnight (asarray)
    f_day = jd + 0.5 - jd05    # Time as fraction of the day (=time/24)
    
    # If jd05 >= jd_start_greg, use the Gregorian calendar:
    sel         = jd05 >= jd_start_greg   # Select the Gregorian cases
    cent_5      = _np.zeros_like(jd05)    # Filled with 0's
    cent_5[sel] = _np.floor((jd05[sel]-1867216.25)/36524.25)   # Count Greg.cent. - 5
    jd05[sel]  += 1 + cent_5[sel] - _np.floor(cent_5[sel]/4.)  # Offset Jul.-Greg.
    
    jd05 += 1524                              # Set JD=0 to -4712-01-01
    yr0   = _np.floor((jd05 - 122.1)/365.25)  # Year since -4716 CE
    day0  = _np.floor(365.25*yr0)             # Julian days since -4716 CE
    mnt0  = _np.floor((jd05-day0)/30.6001)    # Month number + 1 (4-15 = Mar-Feb)
    
    day = jd05 - day0 - _np.floor(30.6001*mnt0) + f_day  # Day of mnt + frac. (as.ar.)
    
    month = _np.zeros_like(mnt0).astype(int)
    month[mnt0 <  14] = (mnt0[mnt0 <  14] -  1)  # Month Mar-Dec: 4-13 -> 3-12
    month[mnt0 >= 14] = (mnt0[mnt0 >= 14] - 13)  # Month Jan,Dec: 14,15 -> 1,2
    
    year = _np.zeros_like(month).astype(int)
    year[month >  2] = (yr0[month >  2] - 4716)  # Year since -4716 CE -> CE
    year[month <= 2] = (yr0[month <= 2] - 4715)
    
    return _np.squeeze(year),_np.squeeze(month),_np.squeeze(day)  # Arrays -> "scalars"



def jd2date_time(jd):
    """Obsolescent function.  Use date_time_from_jd() instead."""
    _warn_obsolescent('jd2date_time', 'date_time_from_jd', rename=True, extra=True)
    return date_time_from_jd(jd)


def date_time_from_jd(jd, jd_start_greg=2299160.5):
    """Compute the calendar date and time from a given Julian Day.
    
    Args:
      jd (float):             Julian day (days).
      jd_start_greg (float):  JD of start of Gregorian calendar (optional; default=2299160.5 = 1582-10-15.0).
    
    Returns:
      tuple (int,int,int, int,int,float):  Tuple containing (year, month, day,  hour, minute, second):
    
        - year (int):      Year CE.  Note that year=0 indicates 1 BCE.
        - month (int):     Month number of year (1-12).
        - day (int):       Day of month (1-31).
    
        - hour (int):      Hour of day (0-23).
        - minute (int):    Minute of hour (0-59).
        - second (float):  Second of minute (0.0-59.999...)
    
    Note:
      - Date and time will be in the same timezone as the JD, hence UTC for the offical JD.
    """
    
    year,month,day_f = date_from_jd(jd, jd_start_greg)  # Day with fraction
    day    = _np.floor(day_f).astype(int)  # Integer day
    time   = (day_f-day)*24                # Time of day in hours
    hour   = _np.floor(time).astype(int)
    minute = _np.floor((time-hour)*60).astype(int)
    second = (time-hour-minute/60)*3600
    
    return year,month,day, hour,minute,second


def jd2year(jd):
    """Obsolescent function.  Use year_from_jd() instead."""
    _warn_obsolescent('jd2year', 'year_from_jd', rename=True)
    return year_from_jd(jd)


def year_from_jd(jd):
    """Compute a year with fraction from a given Julian Day.
    
    Args:
      jd (float):  Julian day (days).
    
    Returns:
      float:  Year CE, with decimals.  Note that year=0 indicates 1 BCE, year=-1 2 BCE, etc.
    """
    
    year,month,day = date_from_jd(jd)         # Compute current year
    ones = _np.ones_like(year)                # Fill "array" with shape of year w. 1s
    jd0  = jd_from_date(year,   ones, ones)   # Jan 1 of current year
    jd1  = jd_from_date(year+1, ones, ones)   # Jan 1 of next year
    dy   = (jd-jd0) / (jd1-jd0)               # Lin. interpol. for fractional year
    
    return year + dy


def fix_date_time(year,month,day, hour,minute,second):
    """Fix a given set of date and time variables (year, month, day, hour, minute and second) to make them
    consistent.
    
    For example, month=13 will be corrected to month=1 in the next year, day=32 to a date in the next month,
    hour=24 to hour=0 on the next day, et cetera.  This is useful, because some sources list hours between 1
    and 24, rather than 0-23, on which Python's datetime crashes.  In rare cases of roundoff of 59.9 or a leap
    second, second=60.  More generally, this makes it straightforward to add or subtract dates and times and
    to take into account timezones, DST, et cetera.
    
    Args:
      year (int):       Year CE.  Note that year=0 = 1 BCE.
      month (int):      Month number of year.
      day (int):        Day of month with fraction.
    
      hour (int):       Hour of time of day.
      minute (int):     Minute of hour of time.
      second (float):   Second of minute of time.
    
    Returns:
      tuple (int,int,int, int,int,float):   tuple containing (year CE, month, day,  hour, minute, second):
    
      - year (int):       Year CE.  Note that year=0 = 1 BCE.
      - month (int):      Month number of year (UT; 1-12).
      - day (int):        Day of month with fraction (UT; 1-31).
    
      - hour (int):       Hour of time of day (0-23).
      - minute (int):     Minute of hour of time (0-59).
      - second (float):   Second of minute of time (0.000-59.999).
        
    Note:
      - uses jd_from_date_time() and date_time_from_jd().
    """
    
    jd = jd_from_date_time(year,month,day, hour,minute,second)
    year,month,day, hour,minute,second = date_time_from_jd(jd)
    
    return year,month,day, hour,minute,second



def doy_from_ymd(year, month, day):
    """Compute the day of year for a given year, month and day.
    
    Args:
      year (int):   Year of date.
      month (int):  Month of date.
      day (int):    Day of date.
    
    Returns:
      (int):        Day of year.
    """
    
    if _np.ndim(year) > 0:  # Array-like:
        # Ensure we have numpy.ndarrays:
        year     = _np.asarray(year)
        month    = _np.asarray(month)
        day      = _np.asarray(day)
        
        ones     = _np.ones(year.shape)  # Array for first month and first day
        
        today    = jd_from_date_time(year, month, _np.floor(day).astype(int),  0, 0, 0.0)
        JanFirst = jd_from_date_time(year, ones, ones,  0, 0, 0.0)
        doy      = _np.floor(today - JanFirst + 1).astype(int)
        
    else:
        today    = jd_from_date_time(year, month, int(_np.floor(day)),  0, 0, 0.0)
        JanFirst = jd_from_date_time(year, 1, 1,  0, 0, 0.0)
        doy      = int(_np.floor(today - JanFirst + 1))
    
    return doy



def doy_from_datetime(date_time):
    """Compute the day of year for a given datetime.
    
    Args:
      date_time (datetime):   Date and time.
      
    Returns:
      (int):        Day of year.
    """
    
    if _np.ndim(date_time) > 0:  # Array-like:
        date_time = _np.asarray(date_time).astype('datetime64[ns]')  # Ensure we have an np.ndarray of type datetime64[ns]
        ymd = ymdhms_us_from_datetime64(date_time)
        doy = doy_from_ymd(ymd[:,0], ymd[:,1], ymd[:,2])
        
    else:  # Scalar:
        doy = doy_from_ymd(date_time.year, date_time.month, date_time.day)
    
    return doy


def month_length_from_year_month(year, month):
    """Return the number of days in a given month.
    
    Parameters:
      year (int):   Year CE (0=1BCE, -1=2BCE, etc.)
      month (int):  Month (1-12).
    
    Returns:
      (int):        Number of days in the desired month.
    """
    
    return date_from_jd( jd_from_date(year,month+1,0))[2]  # Day 0 of next month == last DoM == month length


def jd2tjc(jd):
    """Obsolescent function.  Use tjc_from_jd() instead."""
    _warn_obsolescent('jd2tjc', 'tjc_from_jd', rename=True)
    return tjc_from_jd(jd)


def tjc_from_jd(jd):
    """Compute the time in Julian centuries since 2000.0.
    
    Args:
      jd (float):  Julian day (days).
    
    Returns:
      float:  Time in Julian centuries since 2000.0 (UT).
    """
    
    return (jd - 2451545.0)/36525


def jd2tjm(jd):
    """Obsolescent function.  Use tjm_from_jd() instead."""
    _warn_obsolescent('jd2tjm', 'tjm_from_jd', rename=True)
    return tjm_from_jd(jd)


def tjm_from_jd(jd):
    """Compute the time in Julian millennia since 2000.0.
    
    Args:
      jd (float):  Julian day (days).
    
    Returns:
      float:  Time in Julian millennia since 2000.0 (UT).
    """
    
    return (jd - 2451545.0)/365250


def deltat_1820(jd):
    """Obsolescent function.  Use deltat_from_jd_appr() instead."""
    _warn_obsolescent('deltat_1820', 'deltat_from_jd_appr', rename=True)
    return deltat_from_jd_appr(jd)


def deltat_from_jd_appr(jd):
    """Return a rough approximation for the value of Delta T.
    
    A lenghtening of the day of 1.8 ms/century is assumed, as well as and that the minimum of the parabola is
    DeltaT=12s in 1820.
    
    Args:
      jd (float):  Julian day (days).
    
    Returns:
      float:  Delta T (s).
    
    References:
      - `Extrapolation of Delta T <http://hemel.waarnemen.com/Computing/deltat.html>`_.
    """
    
    # return 12 + 0.5 * 1.8e-3/86400/(36525*86400) * ((jd-jd1820)*86400)**2  # Comprehensible notation
    return 12 + 0.5 * 1.8e-3 / 36525 * (jd-_ac.jd1820)**2                        # Simplified notation


def deltat(jd):
    """Obsolescent function.  Use deltat_from_jd_ipol() instead."""
    _warn_obsolescent('deltat', 'deltat_from_jd_ipol', rename=True)
    return deltat_from_jd_ipol(jd)


def deltat_from_jd_ipol(jd):
    """Return the value of DeltaT through interpolation or 'extrapolation'.
    
    For the date range -700 - now, the value for Delta T is obtained by interpolation of known historical
      values.  Outside this range, a lenghtening of the day of 1.8 ms/century is assumed, as well as that the
      minimum of the parabola is DeltaT=12s in 1820.
    
    Args:
      jd (float):   Julian day (days).
    
    Returns:
      float:  Delta T (s).
    
    References:
      - `International Earth Rotation and Reference Systems Service <ftp://maia.usno.navy.mil/ser7/deltat.data>`_ of the U.S. Naval Observatory.
      - `Robert van Gent's website on Delta T <https://www.staff.science.uu.nl/~gent0113/deltat/deltat.htm>`_.
      - `Extrapolation of Delta T <http://hemel.waarnemen.com/Computing/deltat.html>`_.
    """
    
    # Known data:
    DTyears = [-700,-600,-500,-400,-300,-200,-100,0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1620,1621,1622,1623,1624,1625,1626,1627,1628,1629,1630,1631,1632,1633,1634,1635,1636,1637,1638,1639,1640,1641,1642,1643,1644,1645,1646,1647,1648,1649,1650,1651,1652,1653,1654,1655,1656,1657,1658,1659,1660,1661,1662,1663,1664,1665,1666,1667,1668,1669,1670,1671,1672,1673,1674,1675,1676,1677,1678,1679,1680,1681,1682,1683,1684,1685,1686,1687,1688,1689,1690,1691,1692,1693,1694,1695,1696,1697,1698,1699,1700,1701,1702,1703,1704,1705,1706,1707,1708,1709,1710,1711,1712,1713,1714,1715,1716,1717,1718,1719,1720,1721,1722,1723,1724,1725,1726,1727,1728,1729,1730,1731,1732,1733,1734,1735,1736,1737,1738,1739,1740,1741,1742,1743,1744,1745,1746,1747,1748,1749,1750,1751,1752,1753,1754,1755,1756,1757,1758,1759,1760,1761,1762,1763,1764,1765,1766,1767,1768,1769,1770,1771,1772,1773,1774,1775,1776,1777,1778,1779,1780,1781,1782,1783,1784,1785,1786,1787,1788,1789,1790,1791,1792,1793,1794,1795,1796,1797,1798,1799,1800,1801,1802,1803,1804,1805,1806,1807,1808,1809,1810,1811,1812,1813,1814,1815,1816,1817,1818,1819,1820,1821,1822,1823,1824,1825,1826,1827,1828,1829,1830,1831,1832,1833,1834,1835,1836,1837,1838,1839,1840,1841,1842,1843,1844,1845,1846,1847,1848,1849,1850,1851,1852,1853,1854,1855,1856,1857,1858,1859,1860,1861,1862,1863,1864,1865,1866,1867,1868,1869,1870,1871,1872,1873,1874,1875,1876,1877,1878,1879,1880,1881,1882,1883,1884,1885,1886,1887,1888,1889,1890,1891,1892,1893,1894,1895,1896,1897,1898,1899,1900,1901,1902,1903,1904,1905,1906,1907,1908,1909,1910,1911,1912,1913,1914,1915,1916,1917,1918,1919,1920,1921,1922,1923,1924,1925,1926,1927,1928,1929,1930,1931,1932,1933,1934,1935,1936,1937,1938,1939,1940,1941,1942,1943,1944,1945,1946,1947,1948,1949,1950,1951,1952,1953,1954,1955,1956,1957,1958,1959,1960,1961,1962,1963,1964,1965,1966,1967,1968,1969,1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021]
    DTvalues = [20400,18800,17190,15530,14080,12790,11640,10580,9600,8640,7680,6700,5710,4740,3810,2960,2200,1570,1090,740,490,320,200,120,124,119,115,110,106,102,98,95,91,88,85,82,79,77,74,72,70,67,65,63,62,60,58,57,55,54,53,51,50,49,48,47,46,45,44,43,42,41,40,38,37,36,35,34,33,32,31,30,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,14,13,12,12,11,11,10,10,10,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,14,14,14,14,14,14,14,15,15,15,15,15,15,15,16,16,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,16,16,16,16,15,15,14,14,13.7,13.4,13.1,12.9,12.7,12.6,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.4,12.3,12.2,12.0,11.7,11.4,11.1,10.6,10.2,9.6,9.1,8.6,8.0,7.5,7.0,6.6,6.3,6.0,5.8,5.7,5.6,5.6,5.6,5.7,5.8,5.9,6.1,6.2,6.3,6.5,6.6,6.8,6.9,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.7,7.8,7.8,7.88,7.82,7.54,6.97,6.40,6.02,5.41,4.10,2.92,1.82,1.61,0.10,-1.02,-1.28,-2.69,-3.24,-3.64,-4.54,-4.71,-5.11,-5.40,-5.42,-5.20,-5.46,-5.46,-5.79,-5.63,-5.64,-5.80,-5.66,-5.87,-6.01,-6.19,-6.64,-6.44,-6.47,-6.09,-5.76,-4.66,-3.74,-2.72,-1.54,-0.02,1.24,2.64,3.86,5.37,6.14,7.75,9.13,10.46,11.53,13.36,14.65,16.01,17.20,18.24,19.06,20.25,20.95,21.16,22.25,22.41,23.03,23.49,23.62,23.68,24.49,24.34,24.08,24.02,24.00,23.87,23.95,23.86,23.93,23.73,23.92,23.96,24.02,24.33,24.83,25.30,25.70,26.24,26.77,27.28,27.78,28.25,28.71,29.15,29.57,29.97,30.36,30.72,31.07,31.35,31.68,32.18,32.68,33.15,33.59,34.00,34.47,35.03,35.73,36.54,37.43,38.29,39.20,40.18,41.17,42.23,43.37,44.4841,45.4761,46.4567,47.5214,48.5344,49.5861,50.5387,51.3808,52.1668,52.9565,53.7882,54.3427,54.8712,55.3222,55.8197,56.3000,56.8553,57.5653,58.3092,59.1218,59.9845,60.7853,61.6287,62.2950,62.9659,63.4673,63.8285,64.0908,64.2998,64.4734,64.5736,64.6876,64.8452,65.1464,65.4574,65.7768,66.0699,66.3246,66.6030,66.9069,67.2810,67.6439,68.1024,68.5927,68.9677,69.2202,69.87,70.4]
    
    jd     = _np.asarray(_np.copy(jd))  # Copy and typecast to numpy.ndarray
    jd0    = _np.empty_like(jd)         # Empty array with the shape of jd
    deltat = _np.empty_like(jd)         # Empty array with the shape of jd
    
    year   = year_from_jd(jd)           # Year with fraction for jd of interest
    
    # Cases before the start of our list:
    sel         = year < DTyears[0]
    jd0[sel]    = jd_from_date(DTyears[0],  1, 1)
    deltat[sel] = deltat_from_jd_appr(jd[sel]) - deltat_from_jd_appr(jd0[sel]) + DTvalues[0]
    
    # Cases after the end of our list:
    sel         = year > DTyears[-1]
    jd0[sel]    = jd_from_date(DTyears[-1], 1, 1)
    deltat[sel] = deltat_from_jd_appr(jd[sel]) - deltat_from_jd_appr(jd0[sel]) + DTvalues[-1]
    
    # Cases in our list:
    sel         = _np.logical_and(year >= DTyears[0], year <= DTyears[-1])  # After start, before end
    deltat[sel] = _np.interp(year[sel], DTyears, DTvalues)
    
    return _np.squeeze(deltat)



def gmst(jd):
    """Obsolescent function.  Use gmst_from_jd() instead."""
    _warn_obsolescent('gmst', 'gmst_from_jd', rename=True, extra=True)
    return gmst_from_jd(jd)


def gmst_from_jd(jd, deltat=None):
    """Calculate Greenwich Mean Sidereal Time for any instant, in radians.
    
    Args:
      jd (float):      Julian day (days).
      deltat (float):  Delta T (s).
    
    Returns:
      float:  Greenwich mean sidereal time (rad).
    
    References:
      - Explanatory Supplement to the Astronomical Almanac, 3rd ed, Eq. 6.66 (2012).
    """
    
    tjd  = jd - _ac.jd2000                     # Julian Days after 2000.0 UT
    coefs = [4.89496121042905, 6.30038809894828323, 5.05711849e-15, -4.378e-28, -8.1601415e-29, -2.7445e-36]       # Coefficients for the polynomial
    
    if deltat is None: deltat = 63.8285    # Use DeltaT from J2000 if unavailable
    gmst = 7.0855723730e-12*deltat         # Correction for Delta T
    
    for pow_i,coef_i in enumerate(coefs):  # Add the polynomial
        gmst += coef_i*tjd**pow_i
    
    return gmst % _ac.pi2



def ymdhms_us_from_datetime64(dt64):
    """Convert (array of) datetime64 to a calendar (array of) year, month, day, hour, minute, seconds,
    microsecond with these quantites indexed on the last axis.
    
    Args:
      dt64 (datetime64):  (numpy array of) datetime(s) (of arbitrary shape).
      
    Returns:
       uint32 array:  (..., 7) calendar array with last axis representing year, month, day, hour, minute,
                      second, microsecond.
    
    Note:
      - Nicked from https://stackoverflow.com/a/56260054/1386750
    """
    
    # Allocate output:
    out = _np.empty(dt64.shape + (7,), dtype='u4')
    
    # Decompose calendar floors:
    Y, M, D, h, m, s = [dt64.astype(f'M8[{x}]') for x in 'YMDhms']
    
    out[..., 0] = Y + 1970                     # Gregorian Year
    out[..., 1] = (M - Y) + 1                  # month
    out[..., 2] = (D - M) + 1                  # day
    out[..., 3] = (dt64 - D).astype('m8[h]')   # hour
    out[..., 4] = (dt64 - h).astype('m8[m]')   # minute
    out[..., 5] = (dt64 - m).astype('m8[s]')   # second
    out[..., 6] = (dt64 - s).astype('m8[us]')  # microsecond
    
    return out


def weekday_en_abbr_from_datetime(datetime):
    """Return an English abbreviation of the weekday for a given datetime.
    
    Args:
      datetime (datetime):
    
    Returns:
      (str):  String with the three-character English abbreviation of the weekday (Mon-Sun).
    
    """
    
    weekdays = ['Mon','Tue','Wed','Thu','Fri','Sat','Sun']
    
    return weekdays[datetime.weekday()]


def gps_time_from_jd(jd):
    """Compute GPS time from the Julian day.
    
    Args:
      jd (float):  Julian day (days).
    
    Returns:
      float:  GPS time (seconds since 1980-01-06, without leap seconds!).
    
    Note:
      number of leap seconds on 1980-01-06 is 19.
    """
    
    gps_time = (jd - jd_from_date(1980, 1, 6))*_ac.day_sol + leap_seconds_from_jd(jd) - 19
    
    return gps_time


def jd_from_gps_time(gps_time):
    """Compute the Julian day from the GPS time.
    
    Args:
      gps_time (float):  GPS time (seconds since 1980-01-06, without leap seconds!).
    
    Returns:
      float:  Julian day (days).
    
    Note:
      number of leap seconds on 1980-01-06 is 19.
    """
    
    jd = gps_time/_ac.day_sol + jd_from_date(1980, 1, 6)
    jd -= (leap_seconds_from_jd(jd) - 19)/86400
    
    return jd


def leap_seconds_from_jd(jd, warn=True):
    """Brief description of function.

    Parameters:
      jd (float):   Julian day (days).
      warn (bool):  Warn when offering a date before 1972 (defaults to True).

    Returns:
      (float):  (TAI-UTC), or the number of leap seconds.
    
    Note:
      - Between 1961 and 1972, there were no leap seconds, but continuous functions to compute TAI-UTC.
    
    See:
      - https://hpiers.obspm.fr/eoppc/bul/bulc/UTC-TAI.history
    """
    
    jd        = _np.asarray(_np.copy(jd))  # Copy and typecast to numpy.ndarray
    ls = _np.zeros_like(jd)         # Array with zeros and the shape of jd
    
    if warn & (jd.min() < jd_from_date(1972,1,1)):
        print('Warning: the number of leap seconds is inaccurate before 1972')
    
    mjd = jd - 2400000.5  # Modified Julian day
    
    # Before 1961, there were no leap seconds.  Between 1961 and 1972, there were functions to compute TAI-UTC.
    ls[jd >= jd_from_date(1961,1,1)]  =  1.4228180 + (mjd[jd >= jd_from_date(1961,1,1)]  - 37300) * 0.0012960  # Note that ls~0 around 1958-01-01x
    ls[jd >= jd_from_date(1961,8,1)]  =  1.3728180 + (mjd[jd >= jd_from_date(1961,8,1)]  - 37300) * 0.0012960
    ls[jd >= jd_from_date(1962,1,1)]  =  1.8458580 + (mjd[jd >= jd_from_date(1962,1,1)]  - 37665) * 0.0011232
    ls[jd >= jd_from_date(1963,11,1)] =  1.9458580 + (mjd[jd >= jd_from_date(1963,11,1)] - 37665) * 0.0011232
    ls[jd >= jd_from_date(1964,1,1)]  =  3.2401300 + (mjd[jd >= jd_from_date(1964,1,1)]  - 38761) * 0.0012960
    ls[jd >= jd_from_date(1964,4,1)]  =  3.3401300 + (mjd[jd >= jd_from_date(1964,4,1)]  - 38761) * 0.0012960
    ls[jd >= jd_from_date(1964,9,1)]  =  3.4401300 + (mjd[jd >= jd_from_date(1964,9,1)]  - 38761) * 0.0012960
    ls[jd >= jd_from_date(1965,1,1)]  =  3.5401300 + (mjd[jd >= jd_from_date(1965,1,1)]  - 38761) * 0.0012960
    ls[jd >= jd_from_date(1965,3,1)]  =  3.6401300 + (mjd[jd >= jd_from_date(1965,3,1)]  - 38761) * 0.0012960
    ls[jd >= jd_from_date(1965,7,1)]  =  3.7401300 + (mjd[jd >= jd_from_date(1965,7,1)]  - 38761) * 0.0012960
    ls[jd >= jd_from_date(1965,9,1)]  =  3.8401300 + (mjd[jd >= jd_from_date(1965,9,1)]  - 38761) * 0.0012960
    ls[jd >= jd_from_date(1966,1,1)]  =  4.3131700 + (mjd[jd >= jd_from_date(1966,1,1)]  - 39126) * 0.0025920
    ls[jd >= jd_from_date(1968,2,1)]  =  4.2131700 + (mjd[jd >= jd_from_date(1968,2,1)]  - 39126) * 0.0025920
    
    # Leap seconds:
    ls[jd >= jd_from_date(1972,1,1)] = 10
    ls[jd >= jd_from_date(1972,7,1)] += 1  # Leap second on 1972-07-01
    ls[jd >= jd_from_date(1973,1,1)] += 1  # Leap second on 1973-01-01
    ls[jd >= jd_from_date(1974,1,1)] += 1  # Leap second on 1974-01-01
    ls[jd >= jd_from_date(1975,1,1)] += 1  # Leap second on 1975-01-01
    ls[jd >= jd_from_date(1976,1,1)] += 1  # Leap second on 1976-01-01
    ls[jd >= jd_from_date(1977,1,1)] += 1  # Leap second on 1977-01-01
    ls[jd >= jd_from_date(1978,1,1)] += 1  # Leap second on 1978-01-01
    ls[jd >= jd_from_date(1979,1,1)] += 1  # Leap second on 1979-01-01
    ls[jd >= jd_from_date(1980,1,1)] += 1  # Leap second on 1980-01-01
    ls[jd >= jd_from_date(1981,7,1)] += 1  # Leap second on 1981-07-01
    ls[jd >= jd_from_date(1982,7,1)] += 1  # Leap second on 1982-07-01
    ls[jd >= jd_from_date(1983,7,1)] += 1  # Leap second on 1983-07-01
    ls[jd >= jd_from_date(1985,7,1)] += 1  # Leap second on 1985-07-01
    ls[jd >= jd_from_date(1988,1,1)] += 1  # Leap second on 1988-01-01
    ls[jd >= jd_from_date(1990,1,1)] += 1  # Leap second on 1990-01-01
    ls[jd >= jd_from_date(1991,1,1)] += 1  # Leap second on 1991-01-01
    ls[jd >= jd_from_date(1992,7,1)] += 1  # Leap second on 1992-07-01
    ls[jd >= jd_from_date(1993,7,1)] += 1  # Leap second on 1993-07-01
    ls[jd >= jd_from_date(1994,7,1)] += 1  # Leap second on 1994-07-01
    ls[jd >= jd_from_date(1996,1,1)] += 1  # Leap second on 1996-01-01
    ls[jd >= jd_from_date(1997,7,1)] += 1  # Leap second on 1997-07-01
    ls[jd >= jd_from_date(1999,1,1)] += 1  # Leap second on 1999-01-01
    ls[jd >= jd_from_date(2006,1,1)] += 1  # Leap second on 2006-01-01
    ls[jd >= jd_from_date(2009,1,1)] += 1  # Leap second on 2009-01-01
    ls[jd >= jd_from_date(2012,7,1)] += 1  # Leap second on 2012-07-01
    ls[jd >= jd_from_date(2015,7,1)] += 1  # Leap second on 2015-07-01
    ls[jd >= jd_from_date(2017,1,1)] += 1  # Leap second on 2017-01-01
    
    return _np.squeeze(ls)


def datetime_from_gps_time(gpstime):
    """Create datetime64 objects from GPS times.
    
    Parameters:
      gpstime (float):  Time in GPS-time format.
    
    Returns:
      (numpy.datetime64):  (array of) datetime objects in UTC.
    """
    
    gpstime = _np.asarray(_np.copy(gpstime))          # Copy and typecast to numpy.ndarray
    if gpstime.ndim == 0:  gpstime = gpstime[None]  # Makes gpstime 1D.  Comment: use _np.newaxis instead?
    
    df = _pd.DataFrame()
    df['jd'] = jd_from_gps_time(gpstime)
    df['year'],df['month'],df['day'], df['hour'],df['minute'],df['second'] = date_time_from_jd(df.jd)
    utc = _pd.to_datetime(df[['year','month','day','hour','minute','second']])  # Turn the columns in the df into a single datetime column
    utc = utc.to_numpy()
    
    return _np.squeeze(utc)  # Arrays -> "scalars"


def gps_time_from_datetime(datetimes):
    """Compute GPS times from datetime64 objects.
    
    Parameters:
      datetimes (numpy.datetime64):  (array of) datetime objects in UTC.
    
    Returns:
      (float):  (Array of) time(s) in GPS-time format.
      
    """
    
    datetimes = _np.asarray(_np.copy(datetimes))          # Copy and typecast to numpy.ndarray
    if datetimes.ndim == 0:  datetimes = datetimes[None]  # Makes datetimes 1D.  Comment: use _np.newaxis instead?
    
    jds      = jd_from_datetime(datetimes)
    gpstimes = gps_time_from_jd(jds)
    
    return _np.squeeze(gpstimes)  # Arrays -> "scalars"



def unix_time_from_jd(jd):
    """Compute UNIX time from the Julian day.
    
    Args:
      jd (float):  Julian day (days).
    
    Returns:
      float:  UNIX time (seconds since 1970-01-01).
    """
    
    return (jd - jd_from_date(1970, 1, 1))*_ac.day_sol


def jd_from_unix_time(unix_time):
    """Compute the Julian day from the UNIX time.
    
    Args:
      unix_time (float):  UNIX time (seconds since 1970-01-01).
    
    Returns:
      float:  Julian day (days).
    """
    
    return unix_time/_ac.day_sol + jd_from_date(1970, 1, 1)


def hms_str_from_time(time, use_sec=True, use_ms=False, use_mus=False, use_ns=False):
    """Return a float time in hours as a formatted string in hours, minutes (and seconds).
    
    Parameters:
      time (float):    Time in hours.
      use_sec (bool):  Use seconds in the string.  Optional, defaults to True.
      use_ms (bool):   Use milliseconds in the string.  Optional, defaults to False.
      use_mus (bool):  Use microseconds in the string.  Optional, defaults to False.
      use_ns (bool):   Use nanoseconds in the string.  Optional, defaults to False.
    
    Returns:
      (str):  Time in hours, minutes and seconds of the format hh:mm(:ss(.sss(sss(sss)))).
    """
    
    hr  = _np.floor(time).astype(int)
    mnt = _np.floor((time-hr)*60).astype(int)
    sec = (time-hr-mnt/60)*3600
    
    if use_ns or use_mus or use_ms or (not use_sec):
        if (use_ns and (sec >= 59.9999999995)) or (use_mus and (sec >= 59.9999995)) or \
           (use_ms and (sec >= 59.9995)) or (not use_sec and (sec >= 30)):
            sec = 0
            mnt = mnt+1
    elif (use_sec and (sec >= 59.5)):  # Separate, since use_sec is True by default
        sec = 0
        mnt = mnt+1
        
    if (mnt >= 60):
        mnt = mnt - 60
        hr = hr+1
    if hr >= 24: hr -= 24
    
    if use_ns:
        hms_str = '%2.2i:%2.2i:%012.9f' % (hr, mnt, sec)
    elif use_mus:
        hms_str = '%2.2i:%2.2i:%09.6f' % (hr, mnt, sec)
    elif use_ms:
        hms_str = '%2.2i:%2.2i:%06.3f' % (hr, mnt, sec)
    elif use_sec:
        hms_str = '%2.2i:%2.2i:%2.2i' % (hr, mnt, _np.round(sec,0).astype(int))
    else:
        hms_str = '%2.2i:%2.2i' % (hr, mnt)
        
    return hms_str


def hm_str_from_time(time):
    """Return a float time in hours as a formatted string in hours and minutes: hh:mm.
    
    Parameters:
      time (float):       Time in hours.
    
    Returns:
      (str):  Time in hours and minutes in the format hh:mm.
    
    Note: wrapper for hms_str_from_time().
    """
        
    return hms_str_from_time(time, use_sec=False)


def _warn_obsolescent(old_name, new_name, rename=False, extra=False):
    """Warn that a function is obsolescent and will be removed.  Indicate whether this concerns a simple rename, possibly with extra features."""
    import sys
    sys.stderr.write('\nWarning: the AstroTool function '+old_name+'() is obsolescent and will be removed in a future version.')
    sys.stderr.write('  Use '+new_name+'() instead.')
    if rename:
        if extra:
            sys.stderr.write('  The interface has not changed much; a simple search and replace for the function names should suffice, but please see the documentation for new features.\n\n')
        else:
            sys.stderr.write('  The interface has not changed; a simple search and replace for the function names should suffice.\n\n')
    else:
        sys.stderr.write('  Please see the documentation for details.\n\n')
    return


# Test code:
if __name__ == '__main__':
    import colored_traceback as _clrtrb
    _clrtrb.add_hook()
    
    print('\njd_from_date(), scalar:')
    print('JD for -4712: ', jd_from_date(-4712,1,1.5))
    print('JD for -1000: ', jd_from_date(-1000,1,1.0))
    print('JD for    -1: ', jd_from_date(-1,1,1.0))
    print('JD for     0: ', jd_from_date(0,1,1.0))
    print()
    print('JD for    1J: ', jd_from_date(1,1,1.0))
    print('JD for    1G: ', jd_from_date(1,1,1.0, jd_start_greg=0))
    print()
    print('JD for 1581:  ', jd_from_date(1581,1,1.0))
    print('JD for 1582:  ', jd_from_date(1582,1,1.0))
    print()
    print('JD for 1582a: ', jd_from_date(1582,12,30.0))
    print('JD for 1582b: ', jd_from_date(1582,12,31.0))
    print('JD for 1583a: ', jd_from_date(1583,1,1.0))
    print('JD for 1583b: ', jd_from_date(1583,1,2.0))
    print()
    print('JD for 1584:  ', jd_from_date(1584,1,1.0))
    print('JD for 1585:  ', jd_from_date(1585,1,1.0))
    print()
    print('JD for 2000G: ', jd_from_date(2000,1,1.0))
    print('JD for 2000J: ', jd_from_date(2000,1,1.0, jd_start_greg=_np.inf))
    print()
    print('JD for 2000:  ', jd_from_date_time(2000,1,1.0, 12,0,0))
    
    
    print('\njd_from_date(), array:')
    _years  = [1,1, 1581,1582, 1582,1582,1583,1583, 1584,1585, 2000,2000]
    _months = [1,1,    1,   1,   12,  12,   1,   1,    1,   1,    1,   1]
    _days   = [1,1,    1,   1,   30,  31,   1,   2,    1,   1,    1,   1]
    # julians = [True,False, True,True, True,True,False,False, False,False, False,True]
    _inf = _np.inf
    _jd_start_gregs1 = [_inf,0,  _inf,_inf, _inf,_inf,0,0, 0,0, 0,_inf]
    _letters = ['J:','G:', ': ',': ', ': ',': ',': ',': ', ': ',': ', 'G:','J:']
    
    _JDs = jd_from_date(_years,_months,_days, jd_start_greg=_jd_start_gregs1)
    # _JDs = jd_from_date(_years,_months,_days, _jd_start_greg=_np.inf)
    for _iter in range(len(_JDs)):
        print('%4i%2s  %9.1f' % (_years[_iter],_letters[_iter], _JDs[_iter]))
        
    
    print('\n\ndate_from_jd(), scalar:')
    print('Date for JD=0.0:                    ', *date_from_jd(0.0))
    print('Date for JD=0.0001:                 ', *date_from_jd(0.0001))
    print('Date for JD=1355807.5:              ', *date_from_jd(1355807.5))
    print('Date for JD=1684532.5:              ', *date_from_jd(1684532.5))
    print('Date for JD=1720692.5:              ', *date_from_jd(1720692.5))
    print('Date for JD=1721057.5:              ', *date_from_jd(1721057.5))
    print()
    print('Date for JD=1721423.5:              ', *date_from_jd(1721423.5))
    print('Date for JD=1721423.5 (Gregorian):  ', *date_from_jd(1721423.5, jd_start_greg=0))
    print()
    print('Date for JD=2298152.5:              ', *date_from_jd(2298152.5))
    print('Date for JD=2299160.0:              ', *date_from_jd(2299160.0))
    print('Date for JD=2299161.0:              ', *date_from_jd(2299161.0))
    print('Date for JD=2299238.5:              ', *date_from_jd(2299238.5))
    print('Date for JD=2299247.5:              ', *date_from_jd(2299247.5))
    print('Date for JD=2301795.5:              ', *date_from_jd(2301795.5))
    print()
    print('Date for JD=2451544.5:              ', *date_from_jd(2451544.5))
    print('Date for JD=2451544.5 (Julian):     ', *date_from_jd(2451544.5, jd_start_greg=_np.inf))
    print('Date for JD=2459526.94695:          ', *date_from_jd(2459526.94695))
    print('Date for JD=2459215.54238:          ', *date_from_jd(2459215.54238))
    print()
    
    print('\ndate_from_jd(), array:')
    _jds = _np.arange(26)*1e5
    _jd_start_gregs2 = _np.zeros(len(_jds)) + _np.inf
    # yrs,mnts,dys = date_from_jd(jds)
    _yrs,_mnts,_dys = date_from_jd(_jds, jd_start_greg=_jd_start_gregs2)
    for _iter in range(len(_jds)):
        print('JD = %9.1f:  %5i-%2.2i-%04.1f' % (_jds[_iter], _yrs[_iter], _mnts[_iter], _dys[_iter]))
        
    print()
    for _iter in range(11):
        _jd = 2299160.5-5+_iter
        _yr,_mnt,_dy = date_from_jd(_jd)
        print('JD = %9.1f:  %5i-%2.2i-%04.1f' % (_jd, _yr, _mnt, _dy))
    
    print()
    _yr = 1582
    _mnt = 10
    for _dy in range(11):
        _jd = jd_from_date(_yr,_mnt,_dy)
        print('JD = %9.1f:  %5i-%2.2i-%04.1f' % (_jd, _yr, _mnt, _dy))
    
    
    
    print('\n\ndate_time_from_jd(), scalar:')
    print('Date and time for JD=2459526.94695: ', *date_time_from_jd(2459526.94695))
    print('Date and time for JD=2459215.54238: ', *date_time_from_jd(2459215.54238))
    print()
    print('date_time_from_jd(), array:')
    _jds1 = _jds[15:]
    _yrs,_mnts,_dys, _hrs,_mins,_secs = date_time_from_jd(_jds1+0.127851)  # , jd_start_greg=_jd_start_gregs2)
    for _iter in range(len(_jds1)):
        print('JD = %9.1f:  %5i-%2.2i-%2.2i, %2.2i:%2.2i:%06.3f' % (_jds1[_iter], _yrs[_iter],_mnts[_iter],_dys[_iter], _hrs[_iter],_mins[_iter],_secs[_iter]))
    
    
    
    print('\n\nyear_from_jd(), scalar:')
    print('year for JD=2459526.94695: ', year_from_jd(2459526.94695))
    print('year for JD=2459215.54238: ', year_from_jd(2459215.54238))
    print()
    print('year_from_jd(), array:')
    _jds1 = _jds[15:]
    print(year_from_jd(_jds1+0.127851))
    
    
    print()
    print()
    print('ΔT and GMST for scalars:')
    _jd = jd_from_date_time(2012,12,23, 12,34,56)
    # _jd = jd_from_date_time(2000,1,1, 0,0,0)
    _deltat1 = deltat_from_jd_appr(_jd)
    _deltat2 = deltat_from_jd_ipol(_jd)
    _gmst    = gmst_from_jd(_jd)
    # _gmst1  = gmst_from_jd(_jd, _deltat1)
    # _gmst2  = gmst_from_jd(_jd, _deltat2)
    # print('Date:     ', datetimestr_from_jd(_jd))
    print('JD:       ', _jd)
    print('Delta T1: ', _deltat1,   's')
    print('Delta T2: ', _deltat2,   's')
    print('GMST:     ', _gmst*_ac.r2h,  'h')
    # print('GMST1:    ', _gmst1*_r2h, 'h')
    # print('GMST2:    ', _gmst2*_r2h, 'h')
    
    
    # Leap seconds, GPS times, UNIX times:
    _years  = _np.linspace(1960,2020, 13)  # Every 5 years
    _jds = jd_from_date(_years,1,1)
    
    # GPS times, scalar:
    _gps_time = gps_time_from_jd(_jd)
    _datetime = datetime_from_gps_time(_gps_time)
    print('Leap secs:               ', leap_seconds_from_jd(_jd))
    print('GPS time:                ', _gps_time)
    print('JD from GPS time:        ', jd_from_gps_time(_gps_time))
    print('Datetime:                ', _datetime)
    print('JD from datetime:        ', jd_from_datetime(_datetime))
    print('GPS time from datetime:  ', gps_time_from_datetime(_datetime))
    
    
    # GPS times, arrays:
    print(_jds)
    _gps_times = gps_time_from_jd(_jds)
    _datetimes = datetime_from_gps_time(_gps_times)
    print('Years:                     ', _years)
    print('Leap secs:                 ', leap_seconds_from_jd(_jds))
    print('GPS times:                 ', _gps_times)
    print('JDs from GPS times:        ', jd_from_gps_time(_gps_times))
    print('Datetimes:                 ', _datetimes)
    print('JDs from datetimes:        ', jd_from_datetime(_datetimes))
    print('GPS times from datetimes:  ', gps_time_from_datetime(_datetimes))
    
    _time = 23 + 59/60 + 59.9999999999/3600
    print('hms_str_from_time(): ', hms_str_from_time(_time, use_ns=True))
    print('hms_str_from_time(): ', hms_str_from_time(_time, use_mus=True))
    print('hms_str_from_time(): ', hms_str_from_time(_time, use_ms=True))
    print('hms_str_from_time(): ', hms_str_from_time(_time, use_sec=True))
    print('hms_str_from_time(): ', hms_str_from_time(_time, use_sec=False))
    print('hm_str_from_time():  ', hm_str_from_time(_time))
    
    # UNIX times, scalar:
    _unix_time = unix_time_from_jd(_jd)
    print('UNIX time:   ', _unix_time)
    print('JD:          ', jd_from_unix_time(_unix_time))
    
    # UNIX times, array:
    _unix_times = unix_time_from_jd(_jds)
    print('UNIX times:  ', _unix_times)
    print('JDs:         ', jd_from_unix_time(_unix_times))
