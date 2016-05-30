# A M file to specify the function that generates the doping profile
# With different overlaps and different doping gradient on the drain and source side
# By sebastien Goasguen October 2001/ Purdue

from readinput import *
import numpy as np
import math


def doping(Nx, Ny, Ntotal, junction_l, junction_r, Nd, N_sd, N_body):
    transport_model = transportmodel.value
    fermi_flag = fermiflag1.value

    Lsda = round(Lsd/dx)
    Lg_topa = round(Lg_top/dx)
    Lg_bota = round(Lg_bot/dx)
    t_topa = round(t_top/dy)
    t_bota = round(t_bot/dy)
    t_sia = round(t_si/dy)

    if dopslope_s != 0 or dopslope_d != 0:
        decay_lens = (2/math.sqrt(2.3)*dopslope_s)
        decay_lend = (2/math.sqrt(2.3)*dopslope_d)

    ########################ASSUMING GAUSSIAN DISTRIBUTION########################
    for iii_row in np.arange(((Nx*(t_topa+1))/Nx), ((Ntotal-Nx*t_bota)/Nx-2)):
        for iii_col in np.arange(0, junction_l):
            i_node = iii_row*Nx+iii_col
            Nd[i_node] = N_sd-N_body

        for iii_col in np.arange(junction_l, junction_r-1):
            i_node = iii_row*Nx+iii_col
            if dopslope_s != 0 and dopslope_d != 0:
                Nd[i_node] = N_sd*np.exp(-((iii_col-junction_l+1)*dx/decay_lens)**2)+N_sd*np.exp(-((iii_col-junction_r+1)*dx/decay_lend)**2)-N_body
            elif dopslope_s != 0 and dopslope_d == 0:
                Nd[i_node] = N_sd*np.exp(-((iii_col-junction_l+1)*dx/decay_lens)**2)-N_body
            elif dopslope_s == 0 and dopslope_d != 0:
                Nd[i_node] = N_sd*np.exp(-((iii_col-junction_r)*dx/decay_lend)**2)-N_body

        for iii_col in np.arange(junction_r-1, Nx):
            i_node = iii_row*Nx+iii_col
            Nd[i_node] = N_sd-N_body

    if ox_pnt_flag == 1:
        Nd[(Nx*t_topa):(Nx*t_topa + Nx)] = Nd[(Nx*t_topa + Nx):(Nx*t_topa + 2*Nx)]
        Nd[(Ntotal-Nx*(t_bota+2)):(Ntotal-Nx*(t_bota+2)+Nx)] = Nd[(Ntotal-Nx*(t_bota+2)-Nx):(Ntotal-Nx*(t_bota+2))]

    ######################ABRUPT PROFILE on BOTH SIDE#########################
    elif dopslope_s == 0 and dopslope_d == 0:
        for iii_row in np.arange((Nx*(t_topa+1))/Nx), ((Ntotal-Nx*t_bota)/Nx-2):
            for iii_col in np.arange(0, junction_l):
                i_node = iii_row*Nx+iii_col
                Nd[i_node] = N_sd-N_body

            for iii_col in np.arange(junction_l, junction_r-1):
                i_node = iii_row*Nx+iii_col
                Nd[i_node] = -N_body

            for iii_col in np.arange(junction_r-1, Nx):
                i_node = iii_row*Nx+iii_col
                Nd[i_node] = N_sd-N_body

        if ox_pnt_flag == 1:
            Nd[(Nx*t_topa):(Nx*t_topa) + junction_l] = (N_sd-N_body)/2
            Nd[(Nx*t_topa)+junction_l:(Nx*t_topa-Nx)+junction_r-1] = -N_body/2
            Nd[(Nx*t_topa)+junction_r-1:(Nx*(t_topa+1))] = (N_sd-N_body)/2
            Nd[(Ntotal-Nx*(t_bota+2)):(Ntotal-Nx*(t_bota+2)+junction_l)] = (N_sd-N_body)/2
            Nd[(Ntotal-Nx*(t_bota+2))+junction_l:(Ntotal-Nx*(t_bota+2)+1)+junction_r-1] = -N_body/2
            Nd[(Ntotal-Nx*(t_bota+2))+junction_r-1:(Ntotal-Nx*(t_bota+2)+Nx)] = (N_sd-N_body)/2
    ########################## THE END OF FUNCTION DOPING ######################################################

