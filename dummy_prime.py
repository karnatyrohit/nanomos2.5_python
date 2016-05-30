# Dummy function used to convert quasi Fermi energy to Ne
# (Zhibin Ren 6-5-00)

# See dummy.m for additional comments

from fermi import fermi
import numpy as np


def dummy_prime(x, dummy_flag, fermi_flag):
    if dummy_flag == 0:
        if fermi_flag == 0:
            y = np.exp(x)
        elif fermi_flag == 1:
            y = 1/(1 + np.exp(-x))
    elif dummy_flag == 1/2:
        y = fermi(x, fermi_flag, -1/2)

    return y