from readinput import *
from globvars import globvars
from numpy import linalg
from scipy import sparse
from scipy.sparse.linalg import spsolve
from fermi import fermi

def func_energy(E,tt,U_bias,A,spB_s,spB_d):

    eta = globvars.eta
    Nx = globvars.Nx
    Vd = Vdc.value
    Temp = Te
    fermi_flag = fermiflag1.value

    ee = E
    ep = ee+eta
    ck = 1-((ep-U_bias[0])/(2*tt))
    con_s = -tt*np.exp(1j*np.arccos(ck))
    ck = 1-((ep-U_bias[Nx-1])/(2*tt))
    con_d = -tt*np.exp(1j*np.arccos(ck))
    U_eff = U_bias
    U_eff[0] = U_bias[0]+con_s
    U_eff[Nx-1] = U_bias[Nx-1]+con_d
    G_inv = (ep*np.eye(Nx))-A-np.diag(U_eff)
    G_s = spsolve(sparse.csr_matrix(G_inv), spB_s)
    G_d = spsolve(sparse.csr_matrix(G_inv), spB_d)
    f_1 = fermi(((-Vs-ee)/(k_B*Temp/q)), fermi_flag, -1.0/2.0)
    f_2 = fermi(((-Vd-ee)/(k_B*Temp/q)), fermi_flag, -1.0/2.0)
    N_den = -abs(G_s)**2*np.imag(con_s)*2.0*f_1-abs(G_d)**2*np.imag(con_d)*2.0*f_2
    return N_den