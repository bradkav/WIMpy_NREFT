import numpy as np
import math
from scipy.integrate import simps
from WIMpy import DMUtils as DMU

import pytest


#To-do
#  - Check how the radon transform is normalised!

_Etest = 2.0
_mx = 10.0
_sig = 1e-40

_vlag = 10.0
_sigmav = 156.0
_vesc = 544.0

_Np = 54
_Nn = 131 - _Np


def test_directional_standard():
    theta_list = np.linspace(0, np.pi, 1000)
    dRdtheta_list = 2*np.pi*np.sin(theta_list)*DMU.dRdEdOmega_standard(_Etest, theta_list, _Np, _Nn, _mx, _sig, vlag=_vlag, sigmav=_sigmav, vesc=_vesc)
    integrated_directional = simps(dRdtheta_list, theta_list)
    dRdE = DMU.dRdE_standard(_Etest, _Np, _Nn, _mx, _sig, vlag=_vlag, sigmav=_sigmav, vesc=_vesc)
    assert(math.isclose(integrated_directional, dRdE))
    
@pytest.mark.parametrize('target', DMU.target_list)
def test_directional_NREFT(target):
    cp = 1e-9*np.random.rand(15)
    cn = 1e-9*np.random.rand(15)
    
    theta_list = np.linspace(0, np.pi, 1000)
    dRdtheta_list = 2*np.pi*np.sin(theta_list)*DMU.dRdEdOmega_NREFT(_Etest, theta_list, _mx, cp, cn, target, vlag=_vlag, sigmav=_sigmav, vesc=_vesc)
    integrated_directional = simps(dRdtheta_list, theta_list)
    dRdE = DMU.dRdE_NREFT(_Etest, _mx, cp, cn, target, vlag=_vlag, sigmav=_sigmav, vesc=_vesc)
    
    #print(integrated_directional, dRdE)
    assert(math.isclose(integrated_directional, dRdE))