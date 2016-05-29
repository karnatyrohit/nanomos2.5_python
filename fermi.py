#fermi means integral, the integral order is indicated by the
#last input argument
#x is the input variable, and fermi_flag is the degeneracy flag
#1 for degenerate calculation, 0 for nondegenerate calculation

import numpy as np


def fermi(x, fermi_flag, fermi_order):
    if fermi_order == 1/2:
        if fermi_flag == 1:
            exp_fac = np.exp(-0.17*(x+1.0)**2)
            nu = x**4+50.0+33.6*x*(1.0-0.68*exp_fac)
            zeta=3.0*np.sqrt(np.pi)/(4.0*nu**0.375)
            y=np.exp(x)/(1.0+zeta*np.exp(x))
        elif fermi_flag == 0:
            y=np.exp(x)

    elif fermi_order == 0:
        if fermi_flag == 1:
            y=np.log(1+np.exp(x))
        elif fermi_flag == 0:
            y=np.exp(x)

    elif fermi_order == -1/2:
        if fermi_flag==1:
            exp_fac = np.exp(-0.17*(x+1.0)**2)
            nu = x**4+50.0+33.6*x*(1.0-0.68*exp_fac)
            zeta = 3.0*np.sqrt(np.pi)/(4.0*nu**0.375)
            nu_prime = 4*x**3+33.6-22.848*exp_fac**(1-0.34*(x+x**2))
            zeta_prime = -(9*np.sqrt(np.pi)/32)*nu**(-11/8)*nu_prime
            y = (np.exp(-x)-zeta_prime)/(np.exp(-x)+zeta)**2
        elif fermi_flag == 0:
            y = np.exp(x)

    return y
