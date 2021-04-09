#  Copyright (c) 2019-2021  Marc van der Sluys - marc.vandersluys.nl
#   
#  This file is part of the AstroTool Python package,
#  see: http://astro.ru.nl/~sluys/AstroTool/
#
#  AstroTool has been developed by Marc van der Sluys of the Department of Astrophysics at the Radboud
#  University Nijmegen, the Netherlands and the department of Sustainable energy of the HAN University of
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


"""Date and time functions for AstroTool."""


# Modules:
import numpy as np
from astrotool.constants import pi2, jd1820,jd2000


def julian_day(year,month,day):
    """Compute the Julian Day for a given year, month and day.
    
    Notes:
      - Date and time are expressed in UT.
      - Decimals can be used in the day to take into account the time of day other than midnight, e.g. 1.5 for
        noon on the first day of the month.
    
    Parameters:
      year (int):    Year CE (UT).  Note that year=0 = 1 BCE.
      month (int):   Month number of year (UT; 1-12).
      day (double):  Day of month with fraction (UT; 1.0-31.999).
    
    Returns:
      double:  jd: Julian day (days).
      
    """
    
    if(np.ndim(month) > 0):  # Array-like
        # Ensure we have numpy.ndarrays:
        year  = np.asarray(year)
        month = np.asarray(month)
        day   = np.asarray(day)
        year0 = year
        
        # Jan and Feb are month 13 and 14 of the previous year:
        year[month <= 2] -= 1
        month[month <= 2] += 12
        
        # Assume a Julian date by default:
        a = np.zeros(np.shape(month))
        b = np.zeros(np.shape(month))
        
        # Assume a Gregorian date for years after 1582:
        a[year0 > 1582] = np.floor(year[year0 > 1582]/100.0)
        b[year0 > 1582] = 2 - a[year0 > 1582] + np.floor(a[year0 > 1582]/4.0)
        
    else:                           # If we have a scalar
        
        year0 = year
        
        # Jan and Feb are month 13 and 14 of the previous year:
        if(month <= 2):
            year -= 1
            month += 12
            
            a=0; b=0              # Assume a Julian date by default
            
        if(year0 > 1582):     # Assume a Gregorian date for years after 1582
            a = np.floor(year/100.0)
            b = 2 - a + np.floor(a/4.0)
                
    jd = np.floor(365.25*(year+4716)) + np.floor(30.6001*(month+1)) + day + b - 1524.5
    
    return jd



def date_time2jd(year,month,day, hour,minute,second):
    """Compute the Julian Day for a given year, month, day,  hour, minute and second.
    
    Notes:
      - Date and time are expressed in UT.
    
    Parameters:
      year (int):       Year CE (UT).  Note that year=0 = 1 BCE.
      month (int):      Month number of year (UT; 1-12).
      day (int):        Day of month with fraction (UT; 1.0-31.999).
    
      hour (int):       Hour of time of day (UT).
      minute (int):     Minute of hour of time (UT).
      second (double):  Second of minute of time (UT).
    
    Returns:
      double:  jd: Julian day (days).
    
    Note:
      - uses julian_day().
    """
    
    dday = day + hour/24 + minute/1440 + second/86400
    jd   = julian_day(year, month, dday)
    
    return jd



def jd2ymd(jd):
    """
    Compute the calendar date from a given Julian Day.
    
    Notes:
      - Date and time are expressed in UT.
      - Decimals can be returned in the day to indicate the time of day, e.g. 1.0 for midnight and 1.5 for
        noon on the first day of the month.
    
    Parameters:
      jd (double):  Julian day (days).
    
    Returns:
      tuple (int,int,double):  Tuple containing (year, month, day):
    
        - year (int):    Year CE (UT).  Note that year=0 indicates 1 BCE.
        - month (int):   Month number of year (UT; 1-12).
        - day (double):  Day of month with fraction (UT; 1.0-31.999).
    
    """
    
    z = np.floor(jd+0.5)
    f = jd + 0.5 - z
    
    if(np.ndim(z) > 0):  # Array-like
        z = np.asarray(z)  # Ensure we have a numpy.ndarray
        a = z
        
        # If z >= 2299161, use the Gregorian calendar:
        sel        = (z >= 2299161)
        alpha      = np.zeros(np.shape(z))
        alpha[sel] = np.floor((z[sel]-1867216.25)/36524.25)
        a[sel]     = z[sel] + 1 + alpha[sel] - np.floor(alpha[sel]/4.)
        
    else:                       # If we have a scalar
        
        if(z < 2299161):   # Use the Julian calendar
            a = z
        else:              # Use the Gregorian calendar
            alpha = np.floor((z-1867216.25)/36524.25)
            a = z + 1 + alpha - np.floor(alpha/4.)
    
    
    b = a + 1524
    c = np.floor((b - 122.1)/365.25)
    d = np.floor(365.25*c)
    e = np.floor((b-d)/30.6001)
    day = b - d - np.floor(30.6001*e) + f
    
    if(np.ndim(e) > 0):  # Array-like
        e = np.asarray(e)  # Ensure we have a numpy.ndarray
        month = np.zeros(np.shape(e)).astype(int)
        month[e  < 14] = (e[e  < 14] -  1)
        month[e >= 14] = (e[e >= 14] - 13)
        
        year = np.zeros(np.shape(month)).astype(int)
        year[month  > 2] = (c[month  > 2] - 4716)
        year[month <= 2] = (c[month <= 2] - 4715)
        
    else:                       # If we have a scalar
        if(e < 14):
            month = int(e - 1)
        else:
            month = int(e - 13)
            
        if(month > 2):
            year = int(c - 4716)
        else:
            year = int(c - 4715)
        
    return year,month,day



def jd2year(jd):
    """
    Compute a year with fraction from a given Julian Day.
    
    Parameters:
      jd (double):  Julian day (days).
    
    Returns:
      double:  Year CE, with decimals.  Note that year=0 indicates 1 BCE.
    
    """
    
    year,month,day = jd2ymd(jd)      # Compute current year
    jd0 = julian_day(year,   1, 1)   # Jan 1 of current year
    jd1 = julian_day(year+1, 1, 1)   # Jan 1 of next year
    dy  = (jd-jd0) / (jd1-jd0)       # Linear interpolation for fractional year
    
    return year + dy



def jd2date_time(jd):
    """
    Compute the date and time from a given Julian Day.
    
    Parameters:
      jd (double):  Julian day (days).
    
    Returns:
      tuple (int,int,int, int,int,double):  tuple containing (year CE, month, day,  hour, minute, second):
    
      - year (int):       Year CE (UT).  Note that year=0 = 1 BCE.
      - month (int):      Month number of year (UT; 1-12).
      - day (int):        Day of month with fraction (UT; 1.0-31.999).
    
      - hour (int):       Hour of time of day (UT).
      - minute (int):     Minute of hour of time (UT).
      - second (double):  Second of minute of time (UT).
        
    Note:
      - uses jd2ymd().
    """
    
    year,month,dday = jd2ymd(jd)     # Compute current date
    
    day    = np.floor(dday).astype(int)
    hour   = np.floor((dday-day)*24).astype(int)
    minute = np.floor((dday-day-hour/24)*1440).astype(int)
    second = (dday-day-hour/24-minute/1440)*86400
    
    return year,month,day, hour,minute,second
    


def fix_date_time(year,month,day, hour,minute,second):
    """Fix a given set of date and time variables (year, month, day, hour, minute and second) to make them
    consistent.
    
    For example, month=13 will be corrected to month=1 in the next year, day=32 to a date in the next month,
    hour=24 to hour=0 on the next day, et cetera.  This is useful, because some sources list hours between 1
    and 24, rather than 0-23, on which Python's datetime crashes.  In rare cases of roundoff of 59.9 or a leap
    second, second=60.  More generally, this makes it straightforward to add or subtract dates and times and
    to take into account timezones, DST, et cetera.
    
    Parameters:
      year (int):       Year CE.  Note that year=0 = 1 BCE.
      month (int):      Month number of year.
      day (int):        Day of month with fraction.
    
      hour (int):       Hour of time of day.
      minute (int):     Minute of hour of time.
      second (double):  Second of minute of time.
    
    Returns:
      tuple (int,int,int, int,int,double):  tuple containing (year CE, month, day,  hour, minute, second):
    
      - year (int):       Year CE.  Note that year=0 = 1 BCE.
      - month (int):      Month number of year (UT; 1-12).
      - day (int):        Day of month with fraction (UT; 1-31).
    
      - hour (int):       Hour of time of day (0-23).
      - minute (int):     Minute of hour of time (0-59).
      - second (double):  Second of minute of time (0.000-59.999).
        
    Note:
      - uses date_time2jd() and jd2date_time().
    
    """
    
    jd = date_time2jd(year,month,day, hour,minute,second)
    year,month,day, hour,minute,second = jd2date_time(jd)
    
    return year,month,day, hour,minute,second



def doy_from_ymd(year, month, day):
    """Compute the day of year for a given year, month and day.
    
    Parameters:
      year (int):   Year of date.
      month (int):  Month of date.
      day (int):    Day of date.

    Returns:
      (int):        Day of year.

    """
    
    if(np.ndim(year) > 0):  # Array-like:
        # Ensure we have numpy.ndarrays:
        year     = np.asarray(year)
        month    = np.asarray(month)
        day      = np.asarray(day)
        
        ones     = np.ones(year.shape)  # Array for first month and first day
        
        today    = date_time2jd(year, month, np.floor(day).astype(int),  0, 0, 0.0)
        JanFirst = date_time2jd(year, ones, ones,  0, 0, 0.0)
        doy      = np.floor(today - JanFirst + 1).astype(int)
        
    else:
        today    = date_time2jd(year, month, int(np.floor(day)),  0, 0, 0.0)
        JanFirst = date_time2jd(year, 1, 1,  0, 0, 0.0)
        doy      = int(np.floor(today - JanFirst + 1))
    
    return doy



def doy_from_datetime(date_time):
    """Compute the day of year for a given datetime.
    
    Parameters:
      date_time (datetime):   Date and time.

    Returns:
      (int):        Day of year.

    """
    
    if(np.ndim(date_time) > 0):  # Array-like:
        date_time = np.asarray(date_time)           # Ensure we have an np.ndarray
        ymd = ymdhms_us_from_datetime64(date_time)
        doy = doy_from_ymd(ymd[:,0], ymd[:,1], ymd[:,2])
        
    else:  # Scalar:
        doy = doy_from_ymd(date_time.year, date_time.month, date_time.day)
    
    return doy



def jd2tjc(jd):
    """Compute the time in Julian centuries since 2000.0.
    
    Parameters:
      jd (double):  Julian day (days).
    
    Returns:
      double:  tjc: Time in Julian centuries since 2000.0 (UT).

    """
    
    return (jd - 2451545.0)/36525



def jd2tjm(jd):
    """Compute the time in Julian millennia since 2000.0.
    
    Parameters:
      jd (double):  Julian day (days).
    
    Returns:
      double:  tjm: Time in Julian millennia since 2000.0 (UT).

    """
    
    return (jd - 2451545.0)/365250



def gmst(jd):
    """Calculate Greenwich Mean Sidereal Time for any instant, in radians.
    
    Parameters:
      jd (double):  Julian day (days).
    
    Returns:
      double:  gmst: Greenwich mean sidereal time (rad).
    
    References:
      - Explanatory Supplement to the Astronomical Almanac, 3rd ed, Eq. 6.66 (2012).

    """
    
    tjd  = jd - jd2000                      # Julian Days after 2000.0 UT
    tjd2 = tjd**2
    tjd4 = tjd2**2
           
    gmst = 4.89496121088131 + 6.30038809894828323*tjd + 5.05711849e-15*tjd2 - 4.378e-28*tjd2*tjd - 8.1601415e-29*tjd4 \
        - 2.7445e-36*tjd4*tjd  # Eq. 6.66, removed Delta-T term, hence replaced the first term
    
    return gmst % pi2



def deltat_1820(jd):
    """Return a rough estimate for the value of Delta T.
    
    A lenghtening of the day of 1.8 ms/century is assumed, as well as and that the minimum of the parabola is
    DeltaT=12s in 1820.
    
    Parameters:
      jd (double):  Julian day (days).
    
    Returns:
      double:  Delta T (s).
    
    References:
      - `Extrapolation of Delta T <http://hemel.waarnemen.com/Computing/deltat.html>`_.

    """
    
    # return 12 + 0.5 * 1.8e-3/86400/(36525*86400) * ((jd-jd1820)*86400)**2  # Comprehensible notation
    return 12 + 0.5 * 1.8e-3 / 36525 * (jd-jd1820)**2                        # Simplified notation



def deltat(jd):
    """Return the value of DeltaT through interpolation.
    
    For the date range -700 - now, the value for Delta T is obtained by interpolation of known historical
      values.  Outside this range, a lenghtening of the day of 1.8 ms/century is assumed, as well as that the
      minimum of the parabola is DeltaT=12s in 1820.
    
    Parameters:
      jd (double):  Julian day (days).
    
    Returns:
      double:  Delta T (s).
    
    References:
      - `International Earth Rotation and Reference Systems Service <ftp://maia.usno.navy.mil/ser7/deltat.data>`_ of the U.S. Naval Observatory.
      - `Robert van Gent's website on Delta T <https://www.staff.science.uu.nl/~gent0113/deltat/deltat.htm>`_.
      - `Extrapolation of Delta T <http://hemel.waarnemen.com/Computing/deltat.html>`_.
    
    """
    
    # Known data:
    years = [-700,-600,-500,-400,-300,-200,-100,0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1620,1621,1622,1623,1624,1625,1626,1627,1628,1629,1630,1631,1632,1633,1634,1635,1636,1637,1638,1639,1640,1641,1642,1643,1644,1645,1646,1647,1648,1649,1650,1651,1652,1653,1654,1655,1656,1657,1658,1659,1660,1661,1662,1663,1664,1665,1666,1667,1668,1669,1670,1671,1672,1673,1674,1675,1676,1677,1678,1679,1680,1681,1682,1683,1684,1685,1686,1687,1688,1689,1690,1691,1692,1693,1694,1695,1696,1697,1698,1699,1700,1701,1702,1703,1704,1705,1706,1707,1708,1709,1710,1711,1712,1713,1714,1715,1716,1717,1718,1719,1720,1721,1722,1723,1724,1725,1726,1727,1728,1729,1730,1731,1732,1733,1734,1735,1736,1737,1738,1739,1740,1741,1742,1743,1744,1745,1746,1747,1748,1749,1750,1751,1752,1753,1754,1755,1756,1757,1758,1759,1760,1761,1762,1763,1764,1765,1766,1767,1768,1769,1770,1771,1772,1773,1774,1775,1776,1777,1778,1779,1780,1781,1782,1783,1784,1785,1786,1787,1788,1789,1790,1791,1792,1793,1794,1795,1796,1797,1798,1799,1800,1801,1802,1803,1804,1805,1806,1807,1808,1809,1810,1811,1812,1813,1814,1815,1816,1817,1818,1819,1820,1821,1822,1823,1824,1825,1826,1827,1828,1829,1830,1831,1832,1833,1834,1835,1836,1837,1838,1839,1840,1841,1842,1843,1844,1845,1846,1847,1848,1849,1850,1851,1852,1853,1854,1855,1856,1857,1858,1859,1860,1861,1862,1863,1864,1865,1866,1867,1868,1869,1870,1871,1872,1873,1874,1875,1876,1877,1878,1879,1880,1881,1882,1883,1884,1885,1886,1887,1888,1889,1890,1891,1892,1893,1894,1895,1896,1897,1898,1899,1900,1901,1902,1903,1904,1905,1906,1907,1908,1909,1910,1911,1912,1913,1914,1915,1916,1917,1918,1919,1920,1921,1922,1923,1924,1925,1926,1927,1928,1929,1930,1931,1932,1933,1934,1935,1936,1937,1938,1939,1940,1941,1942,1943,1944,1945,1946,1947,1948,1949,1950,1951,1952,1953,1954,1955,1956,1957,1958,1959,1960,1961,1962,1963,1964,1965,1966,1967,1968,1969,1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021]
    DTvalues = [20400,18800,17190,15530,14080,12790,11640,10580,9600,8640,7680,6700,5710,4740,3810,2960,2200,1570,1090,740,490,320,200,120,124,119,115,110,106,102,98,95,91,88,85,82,79,77,74,72,70,67,65,63,62,60,58,57,55,54,53,51,50,49,48,47,46,45,44,43,42,41,40,38,37,36,35,34,33,32,31,30,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,14,13,12,12,11,11,10,10,10,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,14,14,14,14,14,14,14,15,15,15,15,15,15,15,16,16,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,16,16,16,16,15,15,14,14,13.7,13.4,13.1,12.9,12.7,12.6,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.4,12.3,12.2,12.0,11.7,11.4,11.1,10.6,10.2,9.6,9.1,8.6,8.0,7.5,7.0,6.6,6.3,6.0,5.8,5.7,5.6,5.6,5.6,5.7,5.8,5.9,6.1,6.2,6.3,6.5,6.6,6.8,6.9,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.7,7.8,7.8,7.88,7.82,7.54,6.97,6.40,6.02,5.41,4.10,2.92,1.82,1.61,0.10,-1.02,-1.28,-2.69,-3.24,-3.64,-4.54,-4.71,-5.11,-5.40,-5.42,-5.20,-5.46,-5.46,-5.79,-5.63,-5.64,-5.80,-5.66,-5.87,-6.01,-6.19,-6.64,-6.44,-6.47,-6.09,-5.76,-4.66,-3.74,-2.72,-1.54,-0.02,1.24,2.64,3.86,5.37,6.14,7.75,9.13,10.46,11.53,13.36,14.65,16.01,17.20,18.24,19.06,20.25,20.95,21.16,22.25,22.41,23.03,23.49,23.62,23.68,24.49,24.34,24.08,24.02,24.00,23.87,23.95,23.86,23.93,23.73,23.92,23.96,24.02,24.33,24.83,25.30,25.70,26.24,26.77,27.28,27.78,28.25,28.71,29.15,29.57,29.97,30.36,30.72,31.07,31.35,31.68,32.18,32.68,33.15,33.59,34.00,34.47,35.03,35.73,36.54,37.43,38.29,39.20,40.18,41.17,42.23,43.37,44.4841,45.4761,46.4567,47.5214,48.5344,49.5861,50.5387,51.3808,52.1668,52.9565,53.7882,54.3427,54.8712,55.3222,55.8197,56.3000,56.8553,57.5653,58.3092,59.1218,59.9845,60.7853,61.6287,62.2950,62.9659,63.4673,63.8285,64.0908,64.2998,64.4734,64.5736,64.6876,64.8452,65.1464,65.4574,65.7768,66.0699,66.3246,66.6030,66.9069,67.2810,67.6439,68.1024,68.5927,68.9677,69.2202,69.87,70.4]
    
    year = jd2year(jd)  # Year with fraction for jd of interest
    
    if(year < years[0]):     # Before -700
        jd0 = julian_day(years[0], 1, 1)
        return deltat_1820(jd) - deltat_1820(jd0) + DTvalues[0]
    
    elif(year > years[-1]):  # in the future
        jd1 = julian_day(years[-1], 1, 1)
        return deltat_1820(jd) - deltat_1820(jd1) + DTvalues[-1]
    
    else:                    # linear interpolation from known data
        return np.interp(year, years, DTvalues)



def ymdhms_us_from_datetime64(dt64):
    """Convert (array of) datetime64 to a calendar (array of) year, month, day, hour, minute, seconds,
    microsecond with these quantites indexed on the last axis.
    
    Parameters:
      dt64 (datetime64):  (numpy array of) datetime(s) (of arbitrary shape).
      
    Returns:
       uint32 array:  (..., 7) calendar array with last axis representing year, month, day, hour, minute,
                      second, microsecond.
    
    Note:
      - Nicked from https://stackoverflow.com/a/56260054/1386750
    """
    
    # Allocate output:
    out = np.empty(dt64.shape + (7,), dtype="u4")
    
    # Decompose calendar floors:
    Y, M, D, h, m, s = [dt64.astype(f"M8[{x}]") for x in "YMDhms"]
    
    out[..., 0] = Y + 1970                     # Gregorian Year
    out[..., 1] = (M - Y) + 1                  # month
    out[..., 2] = (D - M) + 1                  # day
    out[..., 3] = (dt64 - D).astype("m8[h]")   # hour
    out[..., 4] = (dt64 - h).astype("m8[m]")   # minute
    out[..., 5] = (dt64 - m).astype("m8[s]")   # second
    out[..., 6] = (dt64 - s).astype("m8[us]")  # microsecond
    
    return out


# Test code:
if(__name__ == "__main__"):
    print(julian_day(2000,1,1.0))
    
