import numpy as np

def Tmin(E_R, A, m_x):
    """
    Minimum DM kinetic energy for recoil of energy E_R
    
    Parameters
    ----------
    * `E_R` [float]:
        Nuclear recoil energy in keV
    * `A` [float]:
        Nuclear mass in atomic mass units
    * `m_x` [float]:
        DM mass in GeV

    
    Returns
    -------
    * `Tmin` [float]:
        Minimum DM kinetic energy in GeV
    
    """
    
    m_N = 0.9315*A
    E_R_GeV = E_R*1e-6
    
    X = np.sqrt(1 + 2*E_R_GeV*(m_x + m_N)**2/(m_N*(2*m_x - E_R_GeV)**2))
    
    sign = np.sign(E_R_GeV - 2*m_x)
    
    Tmin = (E_R_GeV/2 - m_x)*(1 + sign*X)
    
    return Tmin

def vmin_rel(E_R, A, m_x):
    """
    Minimum DM velocity for recoil of energy E_R,
    taking into account relativistic kinematics
    
    Parameters
    ----------
    * `E_R` [float]:
        Nuclear recoil energy in keV
    * `A` [float]:
        Nuclear mass in atomic mass units
    * `m_x` [float]:
        DM mass in GeV
    
    Returns
    -------
    * `vmin_rel` [float]:
        Minimum DM velocity in km/s
    
    """
    
    T = Tmin(E_R, A, m_x)
    
    v = (1 - (T/m_x + 1)**-2)**0.5
    
    return v*3e5