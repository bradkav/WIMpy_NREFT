# coding: utf-8

# DMUtils.py 
# Dark Matter rate calculator as part of WIMpy_NREFT
#
# Author: Bradley J Kavanagh
# Email: bradkav@gmail.com
# Last updated: 24/07/2017

import numpy as np
from numpy import pi, cos, sin
from scipy.integrate import trapz, cumtrapz, quad
from scipy.interpolate import interp1d
from numpy.random import rand
from scipy.special import erf
# from scipy.special import spherical_jn
from Wfunctions import WD, WM, WMP2, WP1, WP2, WS1, WS2, WS1D

#----------------------------------------------------
#---- Velocity Integrals (and helper functions) -----
#----------------------------------------------------


#---------------------------------------------------------
# Velocity integral eta
def calcEta(vmin, vlag=230.0, sigmav=156.0,vesc=544.0):
    aplus = np.minimum((vmin+vlag), vesc)/(np.sqrt(2)*sigmav)
    aminus = np.minimum((vmin-vlag), vesc)/(np.sqrt(2)*sigmav)
    aesc = vesc/(np.sqrt(2)*sigmav)
    
    vel_integral = 0
    
    N = 1.0/(erf(aesc) - np.sqrt(2.0/np.pi)*(vesc/sigmav)*np.exp(-0.5*(vesc/sigmav)**2))
    
    vel_integral = (0.5/vlag)*(erf(aplus) - erf(aminus))
    vel_integral -= (1.0/(np.sqrt(np.pi)*vlag))*(aplus - aminus)*np.exp(-0.5*(vesc/sigmav)**2)
    
    vel_integral = np.clip(vel_integral, 0, 1e30)

    # May need to replace 2pi as below, but probably not
    # return 2*np.pi*N*vel_integral
    return N*vel_integral
    
#---------------------------------------------------------
# Modified velocity integral
def calcMEta(vmin, vlag=230.0, sigmav=156.0,vesc=544.0):
    v0 = np.sqrt(2.0)*sigmav
    amin = vmin/v0
    aplus = np.minimum((vmin+vlag), vesc)/v0
    aminus = np.minimum((vmin-vlag), vesc)/v0
    aesc = vesc/v0
    aE = vlag/v0
    
    A = v0*((aminus/(2*np.sqrt(pi)*aE) + pi**-0.5)*np.exp(-aminus**2) - (aplus/(2*np.sqrt(pi)*aE) - pi**-0.5)*np.exp(-aplus**2))   
    B = (v0/(4.0*aE))*(1+2.0*aE**2)*(erf(aplus) - erf(aminus))
    C = -(v0*pi**-0.5)*(2 + (1/(3.0*aE))*((amin + aesc - aminus)**3 - (amin + aesc - aplus)**3))*np.exp(-aesc**2)
    
    return np.clip(A+B+C, 0, 1e10)/((3e5**2))

#-----------------------------------------------------------
# Minimum velocity 
def vmin(E, m_N, m_x):
    res = E*0.0
    m_N2 = m_N*0.9315
    mu = (m_N2*m_x)/(m_N2+m_x)
    res =  3e5*np.sqrt((E/1e6)*(m_N2)/(2*mu*mu))
    return res
    
#-----------------------------------------------------------
# A helper function for calculating the prefactors to dRdE
def rate_prefactor(m_x):
    rho0 = 0.3
    mu = 1.78e-27*reduced_m(1.0, m_x)
    return 1.38413e-12*rho0/(4.0*np.pi*m_x*mu*mu)
    
    
#-----------------------------------------------------------
# Reduced mass - input A as nucleon number and m_x in GeV
def reduced_m(A, m_x):
    m_A = 0.9315*A
    return (m_A * m_x)/(m_A + m_x)
    
#----------------------------------------------------
#-------------------- Form Factors ------------------
#----------------------------------------------------
    

#-----------------------------------------------------------
# Standard Helm Form Factor for SI scattering
def calcSIFormFactor(E, m_N, old=False):
    #Helm
    if (E < 1e-3):
        return E*0 + 1

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
    return (F**2)*(np.exp(-(q2*s)**2))
    
#--------------------------------------------------------
# Calculate the NREFT nuclear form factors (from 
# coefficients given by 'FFcoeffs'), with couplings cp and cn
# to protons and neutrons
def calcNREFTFormFactor(E,m_A,index,cp,cn, FFcoeffs):
    #Define conversion factor from amu-->keV
    amu = 931.5*1000

    #Convert recoil energy to momentum transfer q in keV
    q1 = np.sqrt(2*m_A*amu*E)
    #Convert q into fm^-1
    q2 = q1*(1e-12/1.97e-7)
    b = np.sqrt(41.467/(45*m_A**(-1.0/3.0) - 25*m_A**(-2.0/3.0)))
    #print " b = ", b, m_A
    y = (q2*b/2)**2

    vpp = FFcoeffs[4*(index-1),:]
    vpn = FFcoeffs[4*(index-1)+1,:]
    vnp = FFcoeffs[4*(index-1)+2,:]
    vnn = FFcoeffs[4*(index-1)+3,:]

    Fpp = np.polyval(vpp[::-1],y)
    Fnn = np.polyval(vnn[::-1],y)
    Fnp = np.polyval(vnp[::-1],y)
    Fpn = np.polyval(vpn[::-1],y)
    return (cn*cn*Fnn + cp*cp*Fpp + cn*cp*Fpn + cp*cn*Fnp)*np.exp(-2*y)

# Replacement Fortran functions
#--------------------------------------------------------
def calcWD(y, target="Xe131", cp=1, cn=1):
    c = np.array([(cp+ cn), (cp- cn)])
    tau1 = np.array([0.,1.])
    tau2 = np.array([0.,1.])
    WDval = []
    for i in range(0,2):
        for j in range(0,2):
            WDval.append(c[i]*c[j]*WD.calcwd(tau1[i], tau2[j], y, target))
    return sum(WDval)

#--------------------------------------------------------
def calcWM(y, target="Xe131", cp=1, cn=1):
    c = np.array([(cp+ cn), (cp- cn)])
    tau1 = np.array([0,1])
    tau2 = np.array([0,1])
    WMval = []
    for i in range(0,2):
        for j in range(0,2):
            WMval.append(c[i]*c[j]*WM.calcwm(tau1[i], tau2[j], y, target))
    return sum(WMval)

#--------------------------------------------------------
def calcWMP2(y, target="Xe131", cp=1, cn=1):
    c = np.array([(cp+ cn), (cp- cn)])
    tau1 = np.array([0.,1.])
    tau2 = np.array([0.,1.])
    WMP2val = []
    for i in range(0,2):
        for j in range(0,2):
            WMP2val.append(c[i]*c[j]*WMP2.calcwmp2(tau1[i], tau2[j], y, target))
    return sum(WMP2val)

#--------------------------------------------------------
def calcWP1(y, target="Xe131", cp=1, cn=1):
    c = np.array([(cp+ cn), (cp- cn)])
    tau1 = np.array([0.,1.])
    tau2 = np.array([0.,1.])
    WP1val = []
    for i in range(0,2):
        for j in range(0,2):
            WP1val.append(c[i]*c[j]*WP1.calcwp1(tau1[i], tau2[j], y, target))
    return sum(WP1val)

#--------------------------------------------------------
def calcWP2(y, target="Xe131", cp=1, cn=1):
    c = np.array([(cp+ cn), (cp- cn)])
    tau1 = np.array([0.,1.])
    tau2 = np.array([0.,1.])
    WP2val = []
    for i in range(0,2):
        for j in range(0,2):
            WP2val.append(c[i]*c[j]*WP2.calcwp2(tau1[i], tau2[j], y, target))
    return sum(WP2val)

#--------------------------------------------------------
def calcWS1(y, target="Xe131", cp=1, cn=1):
    c = np.array([(cp+ cn), (cp- cn)])
    tau1 = np.array([0.,1.])
    tau2 = np.array([0.,1.])
    WS1val = []
    for i in range(0,2):
        for j in range(0,2):
            WS1val.append(c[i]*c[j]*WS1.calcws1(tau1[i], tau2[j], y, target))
    return sum(WS1val)

#--------------------------------------------------------
def calcWS1D(y, target="Xe131", cp=1, cn=1):
    c = np.array([(cp+ cn), (cp- cn)])
    tau1 = np.array([0.,1.])
    tau2 = np.array([0.,1.])
    WS1Dval = []
    for i in range(0,2):
        for j in range(0,2):
            WS1Dval.append(c[i]*c[j]*WS1D.calcws1d(tau1[i], tau2[j], y, target))
    return sum(WS1Dval)

#--------------------------------------------------------
def calcWS2(y, target="Xe131", cp=1, cn=1):
    c = np.array([(cp+ cn), (cp- cn)])
    tau1 = np.array([0.,1.])
    tau2 = np.array([0.,1.])
    WS2val = []
    for i in range(0,2):
        for j in range(0,2):
            WS2val.append(c[i]*c[j]*WS2.calcws2(tau1[i], tau2[j], y, target))
    return sum(WS2val)

# Select the appropriate form factor from the list
#--------------------------------------------------------
def calcFF_M(E, m_A, FFcoeffs, cp=1, cn=1):
  return calcNREFTFormFactor(E, m_A, 1, cp, cn, FFcoeffs)

#--------------------------------------------------------
def calcFF_Sigma1(E, m_A, FFcoeffs, cp=1, cn=1):
  return calcNREFTFormFactor(E, m_A, 2, cp, cn, FFcoeffs)
    
#--------------------------------------------------------
def calcFF_Sigma2(E, m_A,FFcoeffs, cp=1, cn=1):
  return calcNREFTFormFactor(E, m_A, 3, cp, cn, FFcoeffs)   
  
#--------------------------------------------------------
def calcFF_Delta(E, m_A, FFcoeffs,cp=1, cn=1):
  return calcNREFTFormFactor(E, m_A, 4, cp, cn, FFcoeffs)

#--------------------------------------------------------
def calcFF_Phi2(E, m_A, FFcoeffs,cp=1, cn=1):
  return calcNREFTFormFactor(E, m_A, 5, cp, cn, FFcoeffs)

#--------------------------------------------------------
def calcFF_MPhi2(E, m_A, FFcoeffs,cp=1, cn=1):
  return calcNREFTFormFactor(E, m_A, 6, cp, cn, FFcoeffs)

#--------------------------------------------------------
def calcFF_Sigma1Delta(E, m_A, FFcoeffs,cp=1, cn=1):
  return calcNREFTFormFactor(E, m_A, 7, cp, cn, FFcoeffs)



    
#----------------------------------------------
#-------- RECOIL RATES ------------------------
#----------------------------------------------
    
#--------------------------------------------------------
# Standard Spin-Independent recoil rate
# for a particle with (N_p,N_n) protons and neutrons
def dRdE_standard(E, N_p, N_n, mx, sig):
    A = N_p + N_n   
    #print A
    int_factor = sig*calcSIFormFactor(E, A)*(A**2)
    
    return rate_prefactor(mx)*int_factor*calcEta(vmin(E, A, mx), vlag=232, sigmav=156, vesc=544)

#--------------------------------------------------------
# Total number of events for standard SI DM
def Nevents_standard(E_min, E_max, N_p, N_n, mx, sig, eff=None):
    if (eff == None):
        integ = lambda x: dRdE_standard(x, N_p, N_n, mx, sig)
        print(" No efficiency!")
    else:
        integ = lambda x: eff(x)*dRdE_standard(x, N_p, N_n, mx, sig)
    return quad(integ, E_min, E_max)[0]
 
#--------------------------------------------------------
# Differential recoil rate in NREFT framework
# Calculates the contribution from the interference of operators
# i and j (with couplings cp and cn to protons and neutrons)
def dRdE_NREFT_components(E, m_A, m_x, cp, cn, i, j, target, eta, meta, q1):
    #eta = calcEta(vmin(E, m_A, m_x))
    #meta = calcMEta(vmin(E, m_A, m_x))
    amu = 931.5*1000
    #q1 = np.sqrt(2*m_A*amu*E)
    qr = q1/amu
    
    # Required for form factors
    q2 = q1*(1e-12/1.97e-7)
    b = np.sqrt(41.467/(45*m_A**(-1.0/3.0) - 25*m_A**(-2.0/3.0)))
    y = (q2*b/2)**2
        
        
    #Calculate all the form factors, for ease of typing!
    FF_M = lambda x: calcWM(x, target=target, cp=cp, cn=cn)
    FF_Sigma1 = lambda x: calcWS1(x, target=target, cp=cp, cn=cn)
    FF_Sigma2 = lambda x: calcWS2(x, target=target, cp=cp, cn=cn)
    FF_Delta = lambda x: calcWD(x, target=target, cp=cp, cn=cn)
    FF_Phi1 = lambda x: calcWP1(x, target=target, cp=cp, cn=cn)
    FF_Phi2 = lambda x: calcWP2(x, target=target, cp=cp, cn=cn)
    FF_MPhi2 = lambda x: calcWMP2(x, target=target, cp=cp, cn=cn)
    FF_Sigma1Delta = lambda x: calcWS1D(x, target=target, cp=cp, cn=cn)

    #FF_SD = lambda x: calcFF_SD(x, m_A, FFcoeffs, cp, cn)

    rate = 0.0

    #Non-interference terms!
    if (i == j):
        #Contact interactions

        if (i == 1): #STANDARD SPIN-INDEPENDENT
            rate = eta*FF_M(y)
        elif (i == 2):
            rate = 0
        elif (i == 3):
            A = meta*FF_Sigma1(y)
            B = 0.25*(qr**2)*eta*FF_Phi2(y)
            rate = (qr**2)*(A+B)
        elif (i == 4): #STANDARD SPIN-DEPENDENT
            rate = eta*(1.0/16.0)*(FF_Sigma1(y) + FF_Sigma2(y))
            #rate = (1.0/16.0)*eta*FF_SD(E)*(1.0/4.0)
        elif (i == 5):
            A = meta*FF_M(y)
            B = eta*(qr**2)*FF_Delta(y)
            rate = 0.25*(qr**2)*(A+B)
        elif (i == 6):
            rate = (1.0/16.0)*(qr**4)*eta*FF_Sigma2(y)
        elif (i == 7):
            rate =  meta*FF_Sigma1(y)
        elif (i == 8):
            A = meta*FF_M(y)
            B = eta*(qr**2)*FF_Delta(y)
            rate =  0.25*(A+B)
        elif (i == 9):
            rate =  (1.0/16.0)*eta*(qr**2)*FF_Sigma1(y)
        elif (i == 10):
            rate =  0.25*eta*(qr**2)*FF_Sigma2(y)
        elif (i == 11):
            rate =  0.25*eta*(qr**2)*FF_M(y)

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
        



    # We need to do this, because the polynomial form factors
    # aren't valid up to arbitrarily high momenta...
    if (rate < 0):
        return 0.0
    else:
        return rate

#--------------------------------------------------------
# Differential recoil rate in NREFT framework
# Calculates the contribution from the interference of operators
# i and j (with couplings cp and cn to protons and neutrons)
def dRdE_NREFT_sum(E, m_A, m_x, cp, cn, i, j, target, eta, meta, q1):
    #eta = calcEta(vmin(E, m_A, m_x))
    #meta = calcMEta(vmin(E, m_A, m_x))
    amu = 931.5*1000 # keV
    #q1 = np.sqrt(2*m_A*amu*E)
    qr = q1/amu
    
    # Required for form factors
    q2 = q1*(1e-12/1.97e-7)
    b = np.sqrt(41.467/(45*m_A**(-1.0/3.0) - 25*m_A**(-2.0/3.0)))
    y = (q2*b/2)**2
        
        
    #Calculate all the form factors, for ease of typing!
    FF_M = lambda x: calcWM(x, target=target, cp=cp, cn=cn)
    FF_Sigma1 = lambda x: calcWS1(x, target=target, cp=cp, cn=cn)
    FF_Sigma2 = lambda x: calcWS2(x, target=target, cp=cp, cn=cn)
    FF_Delta = lambda x: calcWD(x, target=target, cp=cp, cn=cn)
    FF_Phi1 = lambda x: calcWP1(x, target=target, cp=cp, cn=cn)
    FF_Phi2 = lambda x: calcWP2(x, target=target, cp=cp, cn=cn)
    FF_MPhi2 = lambda x: calcWMP2(x, target=target, cp=cp, cn=cn)
    FF_Sigma1Delta = lambda x: calcWS1D(x, target=target, cp=cp, cn=cn)

    #FF_SD = lambda x: calcFF_SD(x, m_A, FFcoeffs, cp, cn)

    rate = 0.0

    #Non-interference terms!
    if (i == j):
        #Contact interactions

        if (i == 1): #STANDARD SPIN-INDEPENDENT
            rate = eta*FF_M(y)
        elif (i == 2):
            rate = 0
        elif (i == 3):
            A = meta*FF_Sigma1(y)
            B = 0.25*(qr**2)*eta*FF_Phi2(y)
            rate = (qr**2)*(A+B)
        elif (i == 4): #STANDARD SPIN-DEPENDENT
            rate = eta*(1.0/16.0)*(FF_Sigma1(y) + FF_Sigma2(y))
            #rate = (1.0/16.0)*eta*FF_SD(E)*(1.0/4.0)
        elif (i == 5):
            A = meta*FF_M(y)
            B = eta*(qr**2)*FF_Delta(y)
            rate = 0.25*(qr**2)*(A+B)
        elif (i == 6):
            rate = (1.0/16.0)*(qr**4)*eta*FF_Sigma2(y)
        elif (i == 7):
            rate =  meta*FF_Sigma1(y)
        elif (i == 8):
            A = meta*FF_M(y)
            B = eta*(qr**2)*FF_Delta(y)
            rate =  0.25*(A+B)
        elif (i == 9):
            rate =  (1.0/16.0)*eta*(qr**2)*FF_Sigma1(y)
        elif (i == 10):
            rate =  0.25*eta*(qr**2)*FF_Sigma2(y)
        elif (i == 11):
            rate =  0.25*eta*(qr**2)*FF_M(y)

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



    # We need to do this, because the polynomial form factors
    # aren't valid up to arbitrarily high momenta...
    #if (rate < 0):
    #    return 0.0
    #else:
    #    return rate
    return rate
    
#--------------------------------------------------------
# Calculate differential rate in NREFT from two vectors of couplings:
# cp_list and cn_list should be the (dimensional) couplings of the DM
# to protons and neutrons for operators 1-11. That is, cp_list should
# be a vector with 11 elements.
# New version with Form Factors
def dRdE_NREFT(E, m_A, m_x, cp_list, cn_list, target, vlag=230.0, sigmav=156.0,vesc=544.0):
    #Sum over all the contributions from different operators
    
    eta = calcEta(vmin(E, m_A, m_x),vlag=vlag, sigmav=sigmav, vesc=vesc)
    meta = calcMEta(vmin(E, m_A, m_x),vlag=vlag, sigmav=sigmav, vesc=vesc)
    amu = 931.5*1000
    q1 = np.sqrt(2*m_A*amu*E)
    
    dRdE_tot = 0.0
    for i in range(11):
        dRdE_tot += dRdE_NREFT_sum(E, m_A, m_x, cp_list[i], cn_list[i], i+1, i+1, target, eta, meta, q1)
    
    i0 = [1,4,4,8]
    j0 = [3,5,6,9]
    
    for i,j in zip(i0,j0):
        dRdE_tot += 2.0*dRdE_NREFT_sum(E, m_A, m_x, cp_list[i-1], cn_list[j-1], i, j, target, eta, meta, q1)
        
    # conv =  (1.98e-14*1.0/(m_x+amu*1e-6))**2/(16.0*pi)
    conv = (0.3/2./np.pi/m_x)*1.69612985e14 # 1 GeV^-4 * cm^-3 * km^-1 * s * c^6 * hbar^2 to keV^-1 kg^-1 day^-1
    
    # We need to do this, because the polynomial form factors
    # aren't valid up to arbitrarily high momenta...
    dRdE_tot = np.clip(dRdE_tot, 0, 1e30)
    
    return dRdE_tot*conv


#--------------------------------------------------------
# Number of events in NREFT
# See also dRdE_NREFT for more details
# Optionally, you can pass a function 'eff' defining the detector efficiency
def Nevents_NREFT(E_min, E_max, m_A, mx, cp_list, cn_list, target):
    if (eff == None):
        eff = lambda x: 1
    integ = lambda x: eff(x)*dRdE_NREFT(x, m_A, mx, cp_list, cn_list, target)
    
    return quad(integ, E_min, E_max)[0]
    
    
#---------------------------------------------------------
# Code for loading in the FormFactor coefficients
def LoadFormFactors(root, A, Z):
    print(" Loading Form Factor for (A, Z) = (" + str(int(A)) + ", " + str(int(Z)) + ")...")
    return np.loadtxt(root +'/FormFactors_Z='\
                 + str(int(Z)) + '_A=' + str(int(A)) +'.dat')


