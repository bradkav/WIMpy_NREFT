import numpy as np
import math
from scipy.integrate import simps
from WIMpy import DMUtils as DMU


_vlag = 10.0
_sigmav = 156.0
_vesc = 544.0

_vtest = 200.0

#See eq. (I.24) in https://arxiv.org/abs/2104.12785
def test_velocity_integral_normalisation():
    v_list = np.linspace(0, 1000.0, 10000)
    eta0_list = DMU.calcEta(v_list, vlag=_vlag, sigmav=_sigmav,vesc=_vesc)
    assert math.isclose(simps(eta0_list, v_list), 1)
    
#See eq. (I.23) in https://arxiv.org/abs/2104.12785, for n = 0, r = 1
def test_modified_velocity_integral_normalisation():
    v_list = np.linspace(0, 1000.0, 10000)
    eta0_list = DMU.calcEta(v_list, vlag=_vlag, sigmav=_sigmav,vesc=_vesc)
    eta1_norm = DMU.calcMEta(0, vlag=_vlag, sigmav=_sigmav,vesc=_vesc)*(3e5)**2
    assert math.isclose(2*simps(v_list*eta0_list,v_list), eta1_norm)

def test_radon_transform():
    theta_list = np.linspace(0, np.pi, 10000)
    RT_list = DMU.calcRT(_vtest, theta_list, vlag=_vlag, sigmav=_sigmav,vesc=_vesc)
    eta = DMU.calcEta(_vtest, vlag=_vlag, sigmav=_sigmav,vesc=_vesc)
    integrated_RT = simps(np.sin(theta_list)*RT_list, theta_list)
    assert math.isclose(integrated_RT, eta)
    
def test_modified_radon_transform():
    theta_list = np.linspace(0, np.pi, 10000)
    MRT_list = DMU.calcMRT(_vtest, theta_list, vlag=_vlag, sigmav=_sigmav,vesc=_vesc)
    Meta = DMU.calcMEta(_vtest, vlag=_vlag, sigmav=_sigmav,vesc=_vesc)
    integrated_MRT = simps(np.sin(theta_list)*MRT_list, theta_list)
    assert math.isclose(integrated_MRT, Meta)