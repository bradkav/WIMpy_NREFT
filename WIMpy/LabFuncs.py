"""
LabFuncs.py - Version 1 - 27/01/2017

Module for calculating the velocity of a lab
on the Earth as a function of time.

Based on code provided by Ciaran O'Hare.

Author: Bradley J Kavanagh
Email: bradkav@gmail.com
"""

import numpy as np
from numpy import cos, sin, pi, floor

#----Parameter definitions----

vv_earthrev = 29.8
eccentricity = 0.016722
eccentricity_deg = 0.9574
orb_long_ecliptic = 13.0+1.0
v_LSR = 220.0

lat_ecl_gal = np.array([-5.5303,59.575,29.812])
long_ecl_gal = np.array([266.141,-13.3485,179.3212])
v_pec = np.array([11.1,12.2,7.3])


#----Lab Velocity----

#Calculates time in JD for a given date
def JulianDay (month, day, year, hour):
    year_r = year+4800-floor((14-month)/12.0)
    month_r = month+12*floor((14-month)/12.0)-3
    JulianDay = day + floor((153*month_r+2)/5.0) + 365*year_r \
         + floor(year_r/4.0) - floor(year_r/100.0) + \
         floor(year_r/400.0) - 32045 + (hour-12.0)/24.0
    return JulianDay

#Time in Julian days -> v_lab in km/s
def LabVelocity(JD, lat, lon):
    
    #Lab time conversion
    UT = 24*(JD+0.5-floor(JD+0.5)) #Universal time
    MJD = JD - 2400000.5 #Modified Julian Day
    T_0 = (floor(MJD)-55197.5)/36525.0
    t_GAST = (101.0308 + 36000.770*T_0 + 15.04107*UT)/15.0
    t_lab = t_GAST + lon/15
    t_lab = 15*t_lab #Lab time in degrees

    
    #Galactic (LSR) Rotation
    v_galrot1 = np.array([0.0,v_LSR,0.0])
    v_galrot = gal2lab(v_galrot1,t_lab, lat) #transform to lab co-ords
    
    #Peculiar solar Motion
    v_solar1 = v_pec
    v_solar = gal2lab(v_solar1,t_lab, lat) #transform to lab co-ords
    
    #Earth's revolution (first calculate in galactic frame then transform)
    e = eccentricity
    lambda_0 = orb_long_ecliptic
    L = 281.0298 + 36000.77*T_0 + 0.04107*UT
    g = 357.9258 + 35999.05*T_0 + 0.04107*UT
    lambda_sun = L + (1.915 - 0.0048*T_0)*sin(g*pi/180.0) \
         + 0.020*sin(2*g*pi/180.0)
    beta = lat_ecl_gal
    lambda_i = long_ecl_gal
    v_earthrev1 = vv_earthrev*(1-e*sin(pi/180.0*(lambda_sun-lambda_0)))* \
         (cos(beta*pi/180.0)*sin(pi/180.0*(lambda_sun-lambda_i)))
    v_earthrev = gal2lab(v_earthrev1,t_lab, lat) #transform to lab co-ords
    
    #Earth's rotation
    v_earthrot = 0.465102*cos(lat*pi/180)*np.array([0.0,-1.0,0.0]) #already in lab co-ords
    
    #Total
    v_lab = np.array([0.,0.,0.])
    v_lab += v_earthrot
    v_lab += v_earthrev
    v_lab += v_solar
    v_lab += v_galrot
    
    return v_lab
    
    
#----Coordinate transformations----

#Equatorial (x_e,y_e,z_e) to Laboratory (N,W,Z) [TIME DEPENDENT]
def eqt2lab(vp,t_lab, lat):
    t = t_lab*pi/180.0
    latr = lat*pi/180.0
    v = vp*0.0
    v[0] = -cos(t)*sin(latr)*vp[0] - sin(t)*sin(latr)*vp[1] + cos(latr)*vp[2]
    v[1] = sin(t)*vp[0] - cos(t)*vp[1]
    v[2] = cos(t)*cos(latr)*vp[0] + cos(latr)*sin(t)*vp[1] + sin(latr)*vp[2]
    return v

#Galactic (x_g,y_g,z_g) to Equatorial (x_e,y_e,z_e) 
def gal2eqt(vp):
    v = 0.0*vp
    v[0] = -0.06699*vp[0] + 0.4927*vp[1] - 0.8676*vp[2]
    v[1] = -0.8728*vp[0] -0.4503*vp[1] -0.1884*vp[2]
    v[2] = -0.4835*vp[0] + 0.7446*vp[1] + 0.4602*vp[2]
    return v

#Galactic (x_g,y_g,z_g) to Laboratory (N,W,Z) [TIME DEPENDENT]
def gal2lab(v,t_lab, lat):
    vp = gal2eqt(v)
    return eqt2lab(vp, t_lab, lat)