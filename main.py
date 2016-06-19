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
from current import current
from globvars import globvars

def main():
    transport_model = transportmodel.value
    fermi_flag = fermiflag1.value
    Vd = Vdc.value
    Is = globvars.Is

    N_sd = Nsd1.value
    N_body = Nbody1.value

    Lsda = round(Lsd/dx)
    Lg_topa = round(Lg_top/dx)
    Lg_bota = round(Lg_bot/dx)
    t_topa = round(t_top/dy)
    t_bota = round(t_bot/dy)
    t_sia = round(t_si/dy)

    #Parameters for ET model
    delta_T_1 = 5.0/2.0 #energy flux parameter one
    delta_T_2 = 5.0/2.0 #energy flux parameter two
    dim_c=2 #degree of carrier freedom


    ###########################################################################
    # Gate and drain bias layout ##########
    ###########################################################################
    # Calculate total number of bias points
    N_points = (Ng_step+1)*(Nd_step+1)
    print '\nNumber of bias points = %d\n\n' % N_points

    # Gate bias vector
    # Given the number of gate bias steps, step size, and initial gate bias,
    # create a vector containing all gate biases.
    Vg_bias = np.zeros(Ng_step+1)

    # Drain bias vector
    # Given the number of drain bias steps, step size, and initial drain bias,
    # create a vector containing all drain biases.
    Vd_bias = np.zeros(Nd_step+1)


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

    globvars.mx = mx
    globvars.my = my
    globvars.mz = mz

    #########################################################################################
    #SPECIFY THE NEUTRAL BOUNDARY ###########################################################
    #Calculate boundary Ec based neutral charge and Fermi-Dirac statistics###################
    #########################################################################################
    if ox_pnt_flag == 0:
        Nsd1.value = ((t_si/dy)/(t_si/dy-1))*N_sd
        N_sd = Nsd1.value
        Nbody1.value = ((t_si/dy)/(t_si/dy-1))*N_body
        N_body = Nbody1.value

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

    globvars.Nx = Nx

    ###########################################################################
    # Memory allocation
    ###########################################################################
    Ie = np.zeros((Ng_step+1, Nd_step+1))
    Mu_sub_body = np.zeros((t_vall, Ng_step+1, Nd_step+1, Nx, max_subband))
    Ie_sub_body = np.zeros((t_vall, Ng_step+1, Nd_step+1, Nx, max_subband))
    Ne_sub_body = np.zeros((t_vall, Ng_step+1, Nd_step+1, Nx, max_subband))
    Te_sub_body = np.zeros((t_vall, Ng_step+1, Nd_step+1, Nx, max_subband))
    E_sub_body = np.zeros((t_vall, Ng_step+1, Nd_step+1, Nx, max_subband))
    Ne_3d = np.zeros((Nd_step+1, Ntotal, Ng_step+1))
    Ec_3d = np.zeros((Nd_step+1, Ntotal, Ng_step+1))
    conv = {}

    ###############################################################################
    ######################END OF ASSIGNING VARIABLES###############################
    ###############################################################################

    ###############################################################################
    #############################START OF INITIALIZATION###########################
    ###############################################################################
    Nd = np.zeros((Ntotal, 1)) #unchanged through the entire calculation
    F_prime = np.zeros((Ntotal, Ntotal)) #unchanged through the entire
                                   #calculation

    Ne_old = np.zeros((Ntotal, 1))
    Ne_new = np.zeros((Ntotal, 1))
    Ec_old = np.zeros((Ntotal, 1))
    Ec_new = np.zeros((Ntotal, 1))

    Fn_new = np.zeros((Ntotal, 1))

    Ne_sub = np.zeros((t_vall, Nx, max_subband))
    E_sub = np.zeros((t_vall, Nx, max_subband))
    Ne_sub_old = np.zeros((t_vall, Nx, max_subband))
    E_sub_old = np.zeros((t_vall, Nx, max_subband))

    Ie_tem = np.zeros((Nx, 1))
    Ie_sub = np.zeros((t_vall, Nx, max_subband))
    Mu_sub = np.zeros((t_vall, Nx, max_subband))
    Te_sub = np.zeros((t_vall, Nx, max_subband))

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
    nu_scatter = 0
    if transport_model == 5:
        nu_scatter = Nx-2
    elif transport_model == 2:
        nu_scatter = Nx

    globvars.nu_scatter = nu_scatter

    Info_scatter_old = np.zeros((nu_scatter, 4))
    Info_scatter_new = np.zeros((nu_scatter, 4))

    #see reference, MEDICI manual, p2-15
    mu_min = 55*1e-4
    mu_max = 300*1e-4
    Nref = 1e22
    alpha = 0.73

    #============Modified. Mar 18, 2002==================
    Nd2D = np.reshape(Nd,(Ny,Nx)).transpose()
    #============Modified. Mar 18, 2002==================

    for i in np.arange(0, nu_scatter):
        Info_scatter_old[i, 1] = i+1
    #Info_scatter_old(i,4)=1/(1/mu_low+...
    #1/(mu_min+(mu_max-mu_min)./(1+(abs(Nd(Nx*round(t_top/dy)+1+Nx+i))/Nref).^alpha)))
    #Info_scatter_old(i,4)=1/(1/mu_low+1/(mu_min+(mu_max-mu_min)./(1+(abs(Nd2D(i,round(Ny/2)))/Nref).^alpha)))
    #============No Methiessen's rule========================================================
        Info_scatter_old[i, 3] = mu_min+(mu_low-mu_min)/(1+(abs(Nd2D[i, round(Ny/2.0) - 1])/Nref)**alpha)  # Rohit - Check for error due to -1
    #========================================================================================
    globvars.Info_scatter_old = Info_scatter_old
    globvars.Info_scatter_new = Info_scatter_new

    #keyboard
    #############################compress matrix#################################
    spEc = sparse.csr_matrix(Ec_old)
    spNe = sparse.csr_matrix(Ne_old)
    spNd = sparse.lil_matrix(Nd)              #rohit - csr vs lil?
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
    Vdc.value = Vd_initial

    [Fn_new, Ne_new, Ne_sub, E_sub] = charge(spNe, spEc, Ne_sub_old, E_sub_old, Nx, Ny, Ntotal, mx, my, mz, junction_l, junction_r, div_avd)
    #Info_scatter_old=Info_scatter_new
    spFn = sparse.csr_matrix(Fn_new)
    spNe = sparse.csr_matrix(Ne_new)
    Ne_sub_old = Ne_sub
    E_sub_old = E_sub

    [Ec_new] = poisson(Nd, Fn_new, Ec_old, F_prime, div_avd, charge_fac, Eg1, Eg2, Es, Ed, Nx, Ny, Ntotal)  #Rohit - experiment with sparse

    Ec_new = np.reshape(Ec_new, (Ntotal, 1))

    spEc = sparse.csr_matrix(Ec_new)

    transportmodel.value = trans_temp
    print 'transport model = %d' % transportmodel.value
    fermiflag1.value = fermi_temp
    print 'fermi_flag = %d' % fermiflag1.value
    print 'Ntotal = %d' % Ntotal


    if (transport_model!=3) and fermi_flag == 1:
        transport_model = 3
        transportmodel.value = 3
        [Fn_new, Ne_new, Ne_sub, E_sub]=charge(spNe, spEc, Ne_sub_old, E_sub_old, Nx, Ny, Ntotal, mx, my, mz, junction_l, junction_r, div_avd)
        #Info_scatter_old=Info_scatter_new

        spFn = sparse.csr_matrix(Fn_new)
        spNe = sparse.csr_matrix(Ne_new)
        Ne_sub_old = Ne_sub
        E_sub_old = E_sub

        [Ec_new] = poisson(Nd, Fn_new, Ec_new, F_prime, div_avd, charge_fac, Eg1, Eg2, Es, Ed, Nx, Ny, Ntotal)

        Ec_new = np.reshape(Ec_new,(Ntotal,1))
        spEc = sparse.csr_matrix(Ec_new)

        transport_model = trans_temp
        transportmodel.value = trans_temp
        iter_outer = 0
        error_outer = 1
        while error_outer >= criterion_outer:
            [Fn_new,Ne_new,Ne_sub,E_sub] = charge(spNe, spEc, Ne_sub_old, E_sub_old, Nx, Ny, Ntotal, mx, my, mz, junction_l, junction_r, div_avd)
         #Info_scatter_old=Info_scatter_new
            spEc_old = spEc
            spFn = sparse.csr_matrix(Fn_new)
            spNe = sparse.csr_matrix(Ne_new)
            Ne_sub_old = Ne_sub
            E_sub_old = E_sub

            [Ec_new] = poisson(Nd, Fn_new, Ec_new, F_prime, div_avd, charge_fac, Eg1, Eg2, Es, Ed, Nx, Ny, Ntotal)

            Ec_new = np.reshape(Ec_new, (Ntotal, 1))
            spEc = sparse.csr_matrix(Ec_new)
            iter_outer = iter_outer+1
            spEcdiff = spEc - spEc_old
            Ecdiff = spEcdiff.todense()
            error_outer=max(abs((np.real(Ecdiff))))
            print '%s %e \n' % ('error_outer = ', error_outer)

    SpNein = spNe
    SpEcin = spEc
    Ne_sub_oldin = Ne_sub_old
    E_sub_oldin = E_sub_old
    SpNdin = spNd
    SpFnin = spFn
    Fnin = Fn_new
    Ecin = Ec_new
    print 'here'

    ############################END OF INITIAL GUESS OF Ec##############################

    ##########################START OF CURRENT CALCULATION LOOP#########################

    ###############################GATE BIAS LOOP#######################################
    #transport_model=trans_temp;
    Eg1_temp = Eg1
    Eg2_temp = Eg2
    for ii_vg in np.arange (0,Ng_step+1):
        Vg_bias[ii_vg] = Vg1+Vg_step*(ii_vg)
        Eg1 = Eg1_temp - Vg_step*(ii_vg)
        if DG_flag == 1:
            Eg2 = Eg2_temp-Vg_step*(ii_vg)

    # Obtain previous results/initial guess
        spNe=SpNein
        spEc=SpEcin
        Ne_sub_old=Ne_sub_oldin
        E_sub_old=E_sub_oldin
        spNd=SpNdin
        spFn=SpFnin
        Fn_new = Fnin
        Ec_new = Ecin
    ###################################DRAIN BIAS LOOP##################################
        for ii_vd in np.arange(0, Nd_step+1):
            Vd_bias[ii_vd] = Vd_temp+Vd_step*(ii_vd)
            Ed = Ed_temp-Vd_step*(ii_vd)
            Vd = Vd_bias[ii_vd]
            Vdc.value = Vd_bias[ii_vd]
    ############################START OF SELF CONSISTENT LOOP###########################
            iter_outer = 0
            error_outer = 1
            converge = [error_outer]

            while(error_outer>=criterion_outer):

                [Fn_new,Ne_new,Ne_sub,E_sub] = charge(spNe, spEc, Ne_sub_old, E_sub_old, Nx, Ny, Ntotal, mx, my, mz, junction_l, junction_r, div_avd)
                # Info_scatter_old=Info_scatter_new

                spEc_old = spEc
                spFn = sparse.csr_matrix(Fn_new)
                spNe = sparse.csr_matrix(Ne_new)
                Ne_sub_old = Ne_sub
                E_sub_old = E_sub

                [Ec_new] = poisson(Nd,Fn_new, Ec_new, F_prime, div_avd, charge_fac, Eg1, Eg2, Es, Ed, Nx, Ny, Ntotal)  #Rohit - again sparse 

                Ec_new = np.reshape(Ec_new,(Ntotal,1))

                spEc = sparse.csr_matrix(Ec_new)
                iter_outer = iter_outer+1
                spEcdiff2 = spEc - spEc_old
                Ecdiff2 = spEcdiff2.todense()
                error_outer = max(abs(np.real(Ecdiff2)))
                print "iter_outer = %d" % iter_outer
                print 'error_outer = %e' % error_outer

                converge = np.append(converge, [error_outer])

                if iter_outer > 50:
                    ver = '******Converge Problem!!! Please step down DVMAX******'
                    print ver
                    break

            print transport_model
            print 'lalal'
            if transport_model == 5:
                Ie_tem = Is
                #print 'current Is has to be established- coding left'
                #Ie_tem = Is   # Rohit look into this global variable
            else:
                [Ie_tem, Ie_sub, Te_sub, Mu_sub] = current(spNe, spEc, Ne_sub, E_sub, Nx, Ny, Ntotal, mx, my, mz)
    ##########################END OF SELF CONSISTENT LOOP##############################
            Vggg = Vg_bias[ii_vg]
            Vddd = Vd_bias[ii_vd]
            print 'Vggg = %f' % Vggg
            print 'Vddd = %f' % Vddd

            Ie[ii_vg,ii_vd] = np.mean(np.real(Ie_tem))
            Mu_sub_body[:,ii_vg,ii_vd,:,:] = Mu_sub
            Ie_sub_body[:, ii_vg, ii_vd, :, :] = Ie_sub
            Ne_sub_body[:, ii_vg, ii_vd, :, :] = Ne_sub
            Te_sub_body[:, ii_vg, ii_vd, :, :] = Te_sub
            E_sub_body[:, ii_vg, ii_vd, :, :] = E_sub
            Ne_3d[ii_vd, :, ii_vg] = np.reshape(Ne_new, Ntotal)
            Ec_3d[ii_vd, :, ii_vg] = np.reshape(Ec_new, Ntotal)
            conv[ii_vg, ii_vd] = converge

    return [Ie, Ie_sub_body, Te_sub_body, Ne_sub_body, E_sub_body, Ne_3d, Ec_3d, conv, Vd_temp]