# Dummy function used to convert quasi Fermi energy to Ne
# (Zhibin Ren 6-5-00)

# x is the quasi-Fermi energy vector
# dummy_flag indicates whether to use Fermi-Dirac integrals of order zero
#   (2D DOS) or order 1/2 (3D DOS).  This can be related to the thickness of
#   the device and the lattice temperature (see main.m).
# fermi_flag indicates whether to use Fermi-Dirac or Boltzmann statistics

import numpy as np
from fermi import fermi


def dummy(x, dummy_flag, fermi_flag):
    if dummy_flag == 0:
        if fermi_flag == 0:
            y = np.exp(x)
        elif fermi_flag == 1:
            y = np.log(1+np.exp(x))

    elif dummy_flag == 1.0/2.0:
        y = fermi(x, fermi_flag, 1.0/2.0)

    return y