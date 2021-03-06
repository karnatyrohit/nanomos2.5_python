#############################################################################
#######Solving the 1D Schroedinger's equation within vertical slices.########
##########################(Zhibin Ren 7-28-00)###############################
#############################################################################

from readinput import *
from scipy import interpolate
import numpy as np
from numpy import linalg as LA
from scipy import sparse

def schred(Ec_old, Nx, Ny, Ntotal, mx, my, mz):

    transport_model = transportmodel.value
    fermi_flag = fermiflag1.value

    Temp = Te
    N_sd = Nsd1.value
    N_body = Nbody1.value

    #	INPUT AND OUTPUT VARIABLES
    #Ec_old is the old potential energy profile in eV
    #1 column of Ntotal elements
    #E_v: bands formed by subband energies in vertical direction
    #W_v: distribution function in vertical direction

    ############################TEMPORARY VARIABLES##############################

    if ox_pnt_flag==0: # (NO ELECTRON PENETRATION INTO OXIDE REGIONS)
       t_sch = t_si
       Ec_start = Nx*t_top/(dy/refine)
       Ec_end = Ntotal-Nx*t_bot/(dy/refine)
    elif ox_pnt_flag == 1: # (ACCOUNTING FOR ELECTRON PENETRATION INTO OXIDE REGIONS)
       t_sch=t_top+t_si+t_bot
       Ec_start = 0
       Ec_end = Ntotal


    Np_old = round(t_sch/dy)+1
    x_dummy_old = (np.linspace(0, t_sch, Np_old))   # Rohit - verify if transpose is needed?
    Np_new = round(t_sch/(dy/refine))+1
    x_dummy_new = (np.linspace(0, t_sch, Np_new))

    ###############################INITIALIZATION################################

    #print Ec_old
    Ec_old = np.real(Ec_old.todense())
    Ec_old = np.reshape(Ec_old,(Ntotal,1))
    #Ec_old = sparse.csr_matrix(Ec_old)
    E_v = np.zeros((t_vall, Nx, max_subband))
    W_v = np.zeros((max_subband, t_vall, Np_old, Nx))
    W_v_tem_1 = np.zeros((Np_new, 1))
    W_v_tem_2 = np.zeros((Np_old, 1))
    MEc = np.zeros((Np_old,Nx)) # Potential in the silicon region
    Ec_start = round(Ec_start)
    Ec_end = round(Ec_end)
    MEc = (np.reshape(Ec_old[Ec_start:Ec_end], (Np_old, Nx)))
    Ec_old = sparse.csr_matrix(Ec_old)

    if ox_pnt_flag == 0:
        Ec_mod = np.zeros((Np_new,1))
    elif ox_pnt_flag == 1:
        Np_top = round(t_top/(dy/refine))
        Np_bot = round(t_bot/(dy/refine))
        Np_si = round(t_si/(dy/refine))+1
        Ec_top = bar_top*np.ones((Np_top+1, 1))
        Ec_bot = bar_bot*np.ones((Np_bot+1, 1))
        Ec_si = 0*np.ones((Np_si-2, 1))
        Ec_mod = np.array([Ec_top,Ec_si,Ec_bot])

    ##############################################################################
    ################################MAIN COMPUTATION##############################
    ##############################################################################

    for iii_vall in np.arange(0, t_vall):
        m_ee = mz[iii_vall]*m_e
        if iii_vall == 2:
            E_v[2,:,:] = E_v [1,:,:]
            W_v[:, 2, :, :] = W_v[:, 1, :, :]
            break

        tt = (h_bar**2)/(2*m_ee*((dy/refine)**2)*q)

        for iii_col in np.arange(0,Nx):
            if refine == 1.0:
                U_vertical = MEc[:, iii_col]
            else:
                s =interpolate.InterpolatedUnivariateSpline(x_dummy_old, MEc[:,iii_col])
                U_vertical = s(x_dummy_new)
                #U_vertical = interp1(x_dummy_old, MEc[:,iii_col], x_dummy_new, 'spline')

            U_vertical = U_vertical + Ec_mod
            #test = np.diag((U_vertical[1:Np_new-1]).flat)

            H = tt*((2*np.eye(Np_new-2))-(np.diag(np.ones(Np_new-1-2),1))-(np.diag(np.ones(Np_new-1-2),-1))) + np.diag((U_vertical[1:Np_new-1]).flat)
            #print H
            [evalu, evac] = LA.eig(H)
            meval=np.sort(evalu)
            i_order = np.argsort(evalu)
            E_v[iii_vall, iii_col,:] = (meval[0:max_subband])
            for i_counter in np.arange(0,max_subband):
                W_v_tem_1[1:Np_new-1] = np.reshape(np.conjugate(evac[:, i_order[i_counter]]) * evac[:, i_order[i_counter]], (Np_new-2, 1))
            if refine == 1.0:
                W_v_tem_2 = W_v_tem_1
            else:
                s2 = interpolate.InterpolatedUnivariateSpline(x_dummy_new, W_v_tem_1)
                W_v_tem_2 = s2(x_dummy_old)
                # W_v_tem_2 = interp1(x_dummy_new,W_v_tem_1,x_dummy_old,'spline')

            W_v[i_counter, iii_vall, :, iii_col] = np.reshape(W_v_tem_2/sum(W_v_tem_2),Np_new)
    return [E_v, W_v]
    ###########################################################################
    #########################END OF OF SCHRED##################################
    ###########################################################################
