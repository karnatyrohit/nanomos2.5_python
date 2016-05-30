########################################################################################
###################2D_Poisson+1D_Schrodinger+1D_transport##############################
########################################################################################

from readinput import *
import numpy as np
from doping import doping
from fprime import fprime
from charge import charge
from scipy import sparse
from poisson import poisson

def main():
    transport_model = transportmodel.value
    fermi_flag = fermiflag1.value
    global Vd

    N_sd = Nsd1.value
    N_body = Nbody1.value

    Lsda = round(Lsd/dx)
    Lg_topa = round(Lg_top/dx)
    Lg_bota = round(Lg_bot/dx)
    t_topa = round(t_top/dy)
    t_bota = round(t_bot/dy)
    t_sia = round(t_si/dy)

    #Parameters for ET model
    delta_T_1 = 5/2 #energy flux parameter one
    delta_T_2 = 5/2 #energy flux parameter two
    dim_c=2 #degree of carrier freedom

    ##########################################################################################
    ############################Step FUNCTION profile for Nsd#################################
    ##########################################################################################
    junction_l = round((Lsd+overlap_s)/dx)+1
    junction_r = round((Lsd+Lg_top-overlap_d)/dx)+1
    ##########################################################################################
    mx = np.zeros((3, 1)); my = np.zeros((3, 1)); mz = np.zeros((3, 1))
    Temp = Te

    mx[0] = m_t; mx[1] = m_l; mx[2] = m_t
    my[0] = m_t; my[1] = m_t; my[2] = m_l
    mz[0] = m_l; mz[1] = m_t; mz[2] = m_t

    #########################################################################################
    #SPECIFY THE NEUTRAL BOUNDARY ###########################################################
    #Calculate boundary Ec based neutral charge and Fermi-Dirac statistics###################
    #########################################################################################
    if ox_pnt_flag==0:
        Nsd1.value = ((t_si/dy)/(t_si/dy-1))*N_sd
        N_sd = Nsd1.value
        Nbody1.value = ((t_si/dy)/(t_si/dy-1))*N_body
        N_body = Nbody1

    Eg1 = -Vg1+phi_top-psi_si
    Eg2 = -Vg2+phi_bot-psi_si

    if fermi_flag == 0:
        Es = -Vs-k_B*Temp/q*np.log((N_sd-N_body)/Ncc)
        Ed = -Vd-k_B*Temp/q*np.log((N_sd-N_body)/Ncc)
    elif fermi_flag == 1:
        Es = -Vs-k_B*Temp/q*np.log(np.exp((N_sd-N_body)/Ncc)-1)
        Ed = -Vd-k_B*Temp/q*np.log(np.exp((N_sd-N_body)/Ncc)-1)


    #########################################################################################
    ##########################END OF SPECIFY THE NEUTRAL BOUNDARY############################
    #########################################################################################

    ##################################ASSIGNING VARIABLES####################################

    # NORMALIZATION PARAMETERS
    charge_fac = dx*dy*q/(eps_si*eps_o)*Nc
    div_avd = 1e-10*Nc # a parameter used to avoid
    # divergence in converting electron density to dummy quantity

    Nx = round((2*Lsd+Lg_top)/dx)+1
    Ny = round((t_top+t_bot+t_si)/dy)+1
    Ntotal = Nx*Ny

    ###############################################################################
    ######################END OF ASSIGNING VARIABLES###############################
    ###############################################################################

    ###############################################################################
    #############################START OF INITIALIZATION###########################
    ###############################################################################
    Nd = np.zeros(Ntotal, 1) #unchanged through the entire calculation
    F_prime = np.zeros(Ntotal, Ntotal) #unchanged through the entire
                                   #calculation

    Ne_old = np.zeros((Ntotal, 1))
    Ne_new = np.zeros((Ntotal, 1))
    Ec_old = np.zeros((Ntotal, 1))
    Ec_new = np.zeros((Ntotal, 1))

    Fn_new = np.zeros((Ntotal, 1))

    Ne_sub = np.zeros((Nx, max_subband, t_vall))
    E_sub = np.zeros((Nx, max_subband, t_vall))
    Ne_sub_old = np.zeros((Nx, max_subband, t_vall))
    E_sub_old = np.zeros((Nx, max_subband, t_vall))

    Ie_tem = np.zeros((Nx, 1))
    Ie_sub = np.zeros((Nx, max_subband, t_vall))
    Mu_sub = np.zeros((Nx, max_subband, t_vall))
    Te_sub = np.zeros((Nx, max_subband, t_vall))

    ############################START OF SPECIFYING Nd############################
    doping(Nx, Ny, Ntotal, junction_l, junction_r,Nd , N_sd, N_body)
    ###########################END OF SPECIFING Nd###############################

    ###################Preparing F_prime(one time evaluation)####################
    fprime(Nx, Ny, Ntotal, F_prime)
    ###########################END OF SPECIFIING F_prime#########################

    #############################################################################
    #############################END OF INITIALIZATION###########################
    #############################################################################

    #############################################################################
    #############START OF SELF CONSISTENT CALCULATION OF POISSON AND ############
    #############################TRANSPORT EQUATIONS#############################
    #############################################################################
    if transport_model == 5:
        nu_scatter = Nx-2
    elif transport_model == 2:
        nu_scatter = Nx

    Info_scatter_old = np.zeros(nu_scatter,4)
    Info_scatter_new = np.zeros(nu_scatter,4)

    #see reference, MEDICI manual, p2-15
    mu_min = 55*1e-4
    mu_max = 300*1e-4
    Nref = 1e22
    alpha = 0.73

    #============Modified. Mar 18, 2002==================
    Nd2D = np.reshape(Nd,(Ny,Nx)).transpose()
    #============Modified. Mar 18, 2002==================

    for i in range(0, nu_scatter):
        Info_scatter_old[i,2] = i+1
    #Info_scatter_old(i,4)=1/(1/mu_low+...
    #1/(mu_min+(mu_max-mu_min)./(1+(abs(Nd(Nx*round(t_top/dy)+1+Nx+i))/Nref).^alpha)))
    #Info_scatter_old(i,4)=1/(1/mu_low+1/(mu_min+(mu_max-mu_min)./(1+(abs(Nd2D(i,round(Ny/2)))/Nref).^alpha)))
    #============No Methiessen's rule========================================================
        Info_scatter_old[i,4] = mu_min+(mu_low-mu_min)/(1+(abs(Nd2D[i, round(Ny/2) - 1])/Nref)**alpha)  # Rohit - Check for error due to -1
    #========================================================================================

    #keyboard
    #############################compress matrix#################################
    spEc = sparse.csr_matrix(Ec_old)
    spNe = sparse.csr_matrix(Ne_old)
    spNd = sparse.csr_matrix(Nd)
    F_prime = sparse.csr_matrix(F_prime)

    ########################START OF INITIAL GUESS ##############################
    trans_temp = transport_model
    fermi_temp = fermi_flag
    transportmodel.value = 1
    fermiflag1.value = 1
    Vd_temp = Vd
    Ed_temp = Ed
    Ed = Ed_temp+Vd-Vd_initial
    Vd = Vd_initial

    [Fn_new, Ne_new, Ne_sub, E_sub] = charge(spNe, spEc, Ne_sub_old, E_sub_old, Nx, Ny, Ntotal, mx, my, mz, junction_l, junction_r, div_avd)
    #Info_scatter_old=Info_scatter_new
    spFn = sparse.csr_matrix(Fn_new)
    spNe = sparse.csr_matrix(Ne_new)
    Ne_sub_old = Ne_sub
    E_sub_old = E_sub

    [Ec_new] = poisson(spNd,spFn,spEc, F_prime, div_avd, charge_fac, Eg1, Eg2, Es, Ed, Nx, Ny, Ntotal)

    spEc=sparse.csr_matrix(Ec_new)

    transportmodel.value = trans_temp
    fermiflag1.value = fermi_temp
    Ntotal


       if ((transport_model~=3) & fermi_flag==1):

         transport_model=3;

        [Fn_new, Ne_new, Ne_sub, E_sub]=charge(spNe, spEc, Ne_sub_old, E_sub_old, Nx, Ny, Ntotal, mx, my, mz, junction_l, junction_r, div_avd)
       #Info_scatter_old=Info_scatter_new

         spFn=sparse(Fn_new);
         spNe=sparse(Ne_new);
         Ne_sub_old=Ne_sub;
         E_sub_old=E_sub;

         [Ec_new]=poisson(spNd, spFn, spEc, F_prime, div_avd, charge_fac, Eg1, Eg2, Es, Ed, Nx, Ny, Ntotal);

         spEc=sparse(Ec_new);

         transport_model=trans_temp;
         iter_outer=0;
         error_outer=1;
         while(error_outer>=criterion_outer)

          [Fn_new,Ne_new,Ne_sub,E_sub]=charge(spNe,spEc,Ne_sub_old,E_sub_old);
         %Info_scatter_old=Info_scatter_new;

          spEc_old=spEc;
          spFn=sparse(Fn_new);
          spNe=sparse(Ne_new);
          Ne_sub_old=Ne_sub;
          E_sub_old=E_sub;

          [Ec_new]=poisson(spNd, spFn, spEc, F_prime, div_avd, charge_fac, Eg1, Eg2, Es, Ed, Nx, Ny, Ntotal);

          spEc=sparse(Ec_new);
          iter_outer=iter_outer+1;
          error_outer=max(abs(full(real(spEc-spEc_old))));
          fprintf ('%s %e \n','error_outer = ',error_outer);
         end
       end

        SpNein=spNe;
        SpEcin=spEc;
        Ne_sub_oldin=Ne_sub_old;
        E_sub_oldin=E_sub_old;
        SpNdin=spNd;
        SpFnin=spFn;

    ############################END OF INITIAL GUESS OF Ec##############################






    return [Ie,Ie_sub_body,Te_sub_body,Ne_sub_body,E_sub_body,Ne_3d,Ec_3d,conv]