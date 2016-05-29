# anti_dummy function used to convert Ne to quasi Fermi energy
# (Zhibin Ren 6-5-00)

# x is the electron density vector
# dummy_flag
# fermi_flag indicates whether to use Fermi-Dirac or Boltzmann statistics

import numpy as np


def anti_dummy(x, dummy_flag, fermi_flag):
    # Ensure that log(0) won't occur
    xmin = 10e-16
    # replace any value less than xmin with xmin
    x[np.where(x < xmin)] = xmin   # Rohit - verify

    if dummy_flag == 0:
        if fermi_flag == 0:
            y = np.log(x)
        elif fermi_flag == 1:
            y = np.log(np.exp(x)-1)

    elif dummy_flag == 1/2:
        if fermi_flag == 0:
            y = np.log(x)
        elif fermi_flag == 1:
            y = np.log(x) + 3.53553e-1*x - 4.95009e-3*x**2 + 1.48386e-4*x**3 - 4.42563e-6*x**4
    return y
