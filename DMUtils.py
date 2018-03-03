# coding: utf-8

# DMUtils.py 
# Dark Matter rate calculator as part of WIMpy_NREFT
#
# Author: Bradley J Kavanagh
# Email: bradkav@gmail.com
# Last updated: 02/03/2018

import numpy as np
from numpy import pi, cos, sin
from scipy.integrate import trapz, cumtrapz, quad
from scipy.interpolate import interp1d
from numpy.random import rand
from scipy.special import erf

#Nuclear structure functions
from Wfunctions import WD, WM, WMP2, WP1, WP2, WS1, WS2, WS1D


#Load in the list of nuclear spins and atomic masses
target_list = np.loadtxt("Nuclei.txt", usecols=(0,), dtype='string')
J_list = np.loadtxt("Nuclei.txt", usecols=(1,))
A_list = np.loadtxt("Nuclei.txt", usecols=(2,))

#TO-DO:
#   - Implement some kind of 'dict' or something here for the target lists
#   - Replace 'A' as a parameter in the functions, by looking up in A_list


#----------------------------------------------------
#---- Velocity Integrals (and helper functions) -----
#----------------------------------------------------

rho0 = 0.3 #GeV/cm^3

#---------------------------------------------------------
# Velocity integral eta
def calcEta(vmin, vlag=230.0, sigmav=156.0,vesc=544.0):
    
    aplus = np.minimum((vmin+vlag), vmin*0.0 + vesc)/(np.sqrt(2)*sigmav)
    aminus = np.minimum((vmin-vlag), vmin*0.0 + vesc)/(np.sqrt(2)*sigmav)
    aesc = vesc/(np.sqrt(2)*sigmav)
    
    vel_integral = 0
    
    N = 1.0/(erf(aesc) - np.sqrt(2.0/np.pi)*(vesc/sigmav)*np.exp(-0.5*(vesc/sigmav)**2))
    
    vel_integral = (0.5/vlag)*(erf(aplus) - erf(aminus))
    vel_integral -= (1.0/(np.sqrt(np.pi)*vlag))*(aplus - aminus)*np.exp(-0.5*(vesc/sigmav)**2)
    
    vel_integral = np.clip(vel_integral, 0, 1e30)

    return N*vel_integral
    
#---------------------------------------------------------
# Modified velocity integral
def calcMEta(vmin, vlag=230.0, sigmav=156.0,vesc=544.0):
    v0 = np.sqrt(2.0)*sigmav
    amin = vmin/v0
    aplus = np.minimum((vmin+vlag), vmin*0.0 + vesc)/v0
    aminus = np.minimum((vmin-vlag), vmin*0.0 + vesc)/v0
    aesc = vesc/v0
    aE = vlag/v0
    
    A = v0*((aminus/(2*np.sqrt(pi)*aE) + pi**-0.5)*np.exp(-aminus**2) - (aplus/(2*np.sqrt(pi)*aE) - pi**-0.5)*np.exp(-aplus**2))   
    B = (v0/(4.0*aE))*(1+2.0*aE**2)*(erf(aplus) - erf(aminus))
    C = -(v0*pi**-0.5)*(2 + (1/(3.0*aE))*((amin + aesc - aminus)**3 - (amin + aesc - aplus)**3))*np.exp(-aesc**2)
    
    return np.clip(A+B+C, 0, 1e10)/((3e5**2))

#-----------------------------------------------------------
# Minimum velocity 
def vmin(E, A, m_x):
    m_A = A*0.9315
    mu = (m_A*m_x)/(m_A+m_x)
    v =  3e5*np.sqrt((E/1e6)*(m_A)/(2*mu*mu))
    return v
    
#-----------------------------------------------------------
# Reduced mass - input A as nucleon number and m_x in GeV
def reduced_m(A, m_x):
    m_A = 0.9315*A
    return (m_A * m_x)/(m_A + m_x)
    
# A helper function for calculating the prefactors to dRdE
def rate_prefactor(m_x):
    rho0 = 0.3
    mu = 1.78e-27*reduced_m(1.0, m_x)
    return 1.38413e-12*rho0/(m_x*mu*mu)
    
#0.197 GeV  = 1e13/cm
# -> GeV^-1 = 1.97e-14 cm
    
def coupling_to_xsec(c, m_x):
    return (1.97e-14)**2*c**2*(reduced_m(1.0, m_x)**2)/np.pi
    
#----------------------------------------------------
#-------------------- Form Factors ------------------
#----------------------------------------------------
    

#-----------------------------------------------------------
# Standard Helm Form Factor for SI scattering
def calcSIFormFactor(E, m_N, old=False):

    #Define conversion factor from amu-->keV
    amu = 931.5*1e3

    #Convert recoil energy to momentum transfer q in keV
    q1 = np.sqrt(2*m_N*amu*E)

    #Convert q into fm^-1
    q2 = q1*(1e-12/1.97e-7)
    
    #Calculate nuclear parameters
    s = 0.9
    a = 0.52
    c = 1.23*(m_N**(1.0/3.0)) - 0.60
    R1 = np.sqrt(c*c + 7*pi*pi*a*a/3.0 - 5*s*s)
    
    if (old):
        R1 = np.sqrt((1.2**2)*m_N**(2.0/3.0) - 5)
 
    x = q2*R1
    J1 = np.sin(x)/x**2 - np.cos(x)/x
    F = 3*J1/x
    
    formfactor = (F**2)*(np.exp(-(q2*s)**2))
    #formfactor[E < 1e-3] = 1.0
    return formfactor

    
#----------------------------------------------
#-------- RECOIL RATES ------------------------
#----------------------------------------------
    
#--------------------------------------------------------
# Standard Spin-Independent recoil rate
# for a particle with (N_p,N_n) protons and neutrons
def dRdE_standard(E, N_p, N_n, m_x, sig, vlag=232.0, sigmav=156.0, vesc=544.0):
    A = N_p + N_n   
    #print A
    int_factor = sig*calcSIFormFactor(E, A)*(A**2)
    
    return rate_prefactor(m_x)*int_factor*calcEta(vmin(E, A, m_x), vlag, sigmav, vesc)

#--------------------------------------------------------
# Total number of events for standard SI DM
def Nevents_standard(E_min, E_max, N_p, N_n, m_x, sig, eff=None,vlag=232.0, sigmav=156.0, vesc=544.0):
    if (eff == None):
        integ = lambda x: dRdE_standard(x, N_p, N_n, m_x, sig)
        print(" No efficiency!")
    else:
        integ = lambda x: eff(x)*dRdE_standard(x, N_p, N_n, m_x, sig)
    return quad(integ, E_min, E_max)[0]
 
 
#--------------------------------------------------------
# Differential recoil rate in NREFT framework
# Calculates the contribution from the interference of operators
# i and j (with couplings cp and cn to protons and neutrons)
def dRdE_NREFT(E, A, m_x, cp, cn, target, vlag=232.0, sigmav=156.0, vesc=544.0):
    
    eta = calcEta(vmin(E, A, m_x),vlag=vlag, sigmav=sigmav, vesc=vesc)
    meta = calcMEta(vmin(E, A, m_x),vlag=vlag, sigmav=sigmav, vesc=vesc)
    amu = 931.5e3 # keV
    q1 = np.sqrt(2*A*amu*E)

    #Recoil momentum over nucleon mass
    qr = q1/amu
    
    # Required for form factors
    q2 = q1*(1e-12/1.97e-7)
    b = np.sqrt(41.467/(45*A**(-1.0/3.0) - 25*A**(-2.0/3.0)))
    y = (q2*b/2)**2
    
    #Dark matter spin factor
    jx = 0.5
    jfac = jx*(jx+1.0)
    
    rate = E*0.0
    
    c = [cp + cn, cp - cn]
    
    for tau1 in [0,1]:
        for tau2 in [0,1]:
            
            c1 = c[tau1]
            c2 = c[tau2]
    
            R_M = c1[0]*c2[0]*eta + jfac/3.0*(qr**2*meta*c1[4]*c2[4] \
                        + meta*c1[7]*c2[7] + qr**2*eta*c1[10]*c2[10])
            rate += R_M*np.vectorize(WM.calcwm)(tau1, tau2, y, target)
    
            R_P2 = 0.25*qr**2*c1[2]*c2[2]*eta
            rate += R_P2*np.vectorize(WP2.calcwp2)(tau1, tau2, y, target)
    
            #Watch out, this one is the wrong way round...
            R_P2M = eta*c1[2]*c2[0]
            rate += R_P2M*np.vectorize(WMP2.calcwmp2)(tau2, tau1, y, target)
    
            R_S2 = eta*c1[9]*c2[9]*0.25*qr**2 + eta*jfac/12.0*(c1[3]*c2[3] + \
                        qr**2*(c1[3]*c2[5] + c1[5]*c2[3]) + qr**4*c1[5]*c2[5])
            rate += R_S2*np.vectorize(WS2.calcws2)(tau1, tau2, y, target)
    
            R_S1 = (1.0/8.0)*meta*(qr**2*c1[2]*c2[2] + c1[6]*c2[6]) +\
                        jfac/12.0*eta*(c1[3]*c2[3] + qr**2*c1[8]*c2[8])
            rate += R_S1*np.vectorize(WS1.calcws1)(tau1, tau2, y, target)
    
            R_D = jfac/3.0*eta*(qr**2*c1[4]*c2[4] + c1[7]*c2[7])
            rate += R_D*np.vectorize(WD.calcwd)(tau1, tau2, y, target)
    
            #This one might be flipped too
            R_S1D = jfac/3.0*eta*(c1[4]*c2[3] - c1[7]*c2[8])
            rate += R_S1D*np.vectorize(WS1D.calcws1d)(tau2, tau1, y, target)

    conv = (rho0/2./np.pi/m_x)*1.69612985e14 # 1 GeV^-4 * cm^-3 * km^-1 * s * c^6 * hbar^2 to keV^-1 kg^-1 day^-1

    rate = np.clip(rate, 0, 1e30)
    return (4*np.pi/(2*(J_list[target_list == target])+1))*rate*conv


#--------------------------------------------------------
# Number of events in NREFT
# See also dRdE_NREFT for more details
# Optionally, you can pass a function 'eff' defining the detector efficiency
def Nevents_NREFT(E_min, E_max, A, m_x, cp, cn, target, eff = None,vlag=232.0, sigmav=156.0, vesc=544.0):
    if (eff == None):
        eff = lambda x: 1
    integ = lambda x: eff(x)*dRdE_NREFT(x, A, m_x, cp, cn, target, vlag, sigmav, vesc)
    
    return quad(integ, E_min, E_max)[0]



#---------------------------------------------------------
#Code for the long range interactions (which we're not using...)
"""
        #Long-range interactions
        elif (i == 101):
            rate =  (qr**-4)*eta*FF_M(y)
        elif (i == 104):
            #rate =  (qr**-4)*(1.0/16.0)*eta*FF_SD(E)
            rate = 0    #ZERO BY DEFINITION!
        elif (i == 105):
            A = meta*FF_M(y)
            B = eta*(qr**2)*FF_Delta(y)
            rate =  0.25*(qr**-2.0)*(A+B)
        elif (i == 106):
            rate =  (1.0/16.0)*eta*FF_Sigma2(y)
        elif (i == 111):
            rate =  0.25*eta*(qr**-2)*FF_M(y)


    #Interference terms
    else:
        if ((i == 1 and j == 3) or (i == 3 and j == 1)):
            rate = (1.0/2.0)*(qr**2)*eta*FF_MPhi2(y)
        elif ((i == 4 and j == 5) or (i == 5 and j == 4)):
            rate = -(1.0/8.0)*(qr**2)*eta*FF_Sigma1Delta(y)
        elif ((i == 4 and j == 6) or (i == 6 and j == 4)):
            rate = (1.0/16.0)*(qr**2)*eta*FF_Sigma2(y)
        elif ((i == 8 and j == 9) or (i == 9 and j ==8)):
            rate =  (1.0/8.0)*(qr**2)*eta*FF_Sigma1Delta(y)
        elif ((i == 104 and j == 105) or (i == 105 and j == 104)):
            rate =  -(1.0/8.0)*eta*FF_Sigma1Delta(y)
        elif ((i == 104) and (j == 106) or (i == 106 and j == 104)):
            rate =  (1.0/16.0)*eta*FF_Sigma2(y)

"""
