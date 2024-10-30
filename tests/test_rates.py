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


@pytest.mark.parametrize('target', DMU.target_list)
def test_spin_independent(target):
    cp = np.zeros(15)
    cp[0] = 1e-9*np.random.rand(1)
    cn = 1.0*cp
    
    sig = DMU.coupling_to_xsec(cp[0], _mx)
    A = DMU.Avals[target]
    
    E_list = np.geomspace(1e-3, 100)
    
    dRdE_SI = DMU.dRdE_standard(E_list, A, 0 , _mx, sig)
    dRdE_NR = DMU.dRdE_NREFT(E_list, _mx, cp, cn, target)
    
    #mask = dRdE_SI > 0
    
    #errors = ((dRdE_SI[mask] - dRdE_NR[mask])/dRdE_SI[mask])**2
    np.testing.assert_almost_equal(dRdE_SI, dRdE_NR, decimal=6)
    



    