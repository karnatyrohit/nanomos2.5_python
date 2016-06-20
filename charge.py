##############################################################################
###########1D transport equation slover(Zhibin Ren 8-02-00)###################
##############################################################################

from readinput import *
import numpy as np
from schred import schred
from fermi import fermi
from integral import integral
from anti_dummy import anti_dummy
from scipy import sparse
from myquad import myquad
from func_energy import func_energy
from globvars import globvars
from scipy.sparse.linalg import spsolve
from current_mat import current_mat
import cmath


def charge(Ne_old,Ec_old,Ne_sub_old,E_sub_old, Nx, Ny, Ntotal, mx, my, mz, junction_l, junction_r, div_avd):

    transport_model = transportmodel.value
    fermi_flag = fermiflag1.value
    Vd = Vdc.value
    N_dos = globvars.N_dos

    Lsda=round(Lsd/dx)
    Lg_topa=round(Lg_top/dx)
    Lg_bota=round(Lg_bot/dx)
    t_topa=round(t_top/dy)
    t_bota=round(t_bot/dy)
    t_sia=round(t_si/dy)
    Temp = Te

    #NEGF scattering method parameters

    nu_scatter = globvars.nu_scatter
    Info_scatter_new  = globvars.Info_scatter_new
    Info_scatter_old = globvars.Info_scatter_old
    Is = globvars.Is
    Id = globvars.Id

    #######################MODEL 4 related parameters#############################
    eta = 1e-6j
    globvars.eta = 1e-6j
    Ef_tail_low = 0.01
    Ef_tail_up = 0.3
    E_step = criterion_outer/2.0 #0.5 times criterion_outer
    zeta_self = -(25.0*120.0/(mu_low*1e4))*1.0e-3j #eV
    decay_fac = 1 #dimensionless
    ########################MODEL 3 related parameters############################
    ########################MODEL 2 related parameters############################
    e_field_limit = 1.0 #in V/m
    ##############################################################################
    ###########################INPUT AND OUTPUT VARIABLES#########################
    ##############################################################################
    # Ne_old is the old electron density
    # 1 column of Ntotal elements
    # Ec_old is the old potential energy profile in eV
    # 1 column of Ntotal elements
    # Te_old is the old electron temperature along the channel
    # 1 column of Nx elements, for charge-sheet-model

    # Ne_sub_old contains 2D electron density on each subband
    # E_sub_old contains the 1D subband energy for all subband
    # considered
    # Te_sub_old contains temperture data of electron on each subband
    # Nx by m, for bulk model

    # Fn_new 1 column of Ntotal elements, dummy varible
    # for Poisson Equation solving
    # Ne_new is the new 3D density in m^-3
    # 1 column of Ntotal elements
    # Ne_sub contains 2D electron density on each subband
    # E_sub contains the 1D subband energy for all subband considered

    # Fn is dummy Fermi energy level used by models 1,4 and 5 for
    # current calculation
    # max_subband is the number of subbands considered
    # t_vall is the number of ladders considered
    #############################################################################

    #########################INITIALIZATION######################################
    Ec_old = np.real(Ec_old.todense())#1 column of Ntotal elements
    Ec_old = sparse.csr_matrix(Ec_old)
    Ne_old = np.real(Ne_old.todense()) #1 column of Ntotal elements
    #Ne_old = sparse.csr_matrix(Ne_old)
    #Ne_sub_old = zeros(Nx,max_subband,t_vall)
    #E_sub_old = zeros(Nx,max_subband,t_vall)

    Fn_new = np.zeros((Ntotal, 1)) #1 column of Ntotal elements
    Ne_new = np.zeros((Ntotal, 1)) #1 column of Ntotal elements
    Ne_sub = np.zeros((t_vall, Nx, max_subband))
    E_sub = np.zeros((t_vall, Nx, max_subband))

    # MODIFIED BY RAMESH,  JAN 13th 03. Comment out as
    # already initialized in main.m
    # Info_scatter_new = np.zeros(nu_scatter, 4)
    # Info_scatter_new = Info_scatter_old

    # temporarily used variables
    if ox_pnt_flag == 0:
        Np_v = round(t_si/dy)+1
    elif ox_pnt_flag == 1:
        Np_v = Ny


    Fn = np.zeros((Nx, 1))
    Ns = np.zeros((Nx, 1))
    N_body = np.zeros((Np_v, Nx))
    N_body_sum = np.zeros((Np_v, Nx))
    U_sub = np.zeros((t_vall, Nx, max_subband))
    W_sub = np.zeros((max_subband, t_vall, Np_v, Nx))
    U_bias = np.zeros((Nx, 1))

    ############################################################################
    #THE START OF TRANSPORT EQUATION PART#######################################
    ############################################################################

    ############################################################################
    if transport_model == 1:  # TRANSPORT MODEL 1####################################
    ############################################################################
    # INITIALIZE THE REAL FERMI ENERGY LEVEL AND COMPUTE
        N_row_mid = round((Ntotal-Nx*(t_topa+t_bota+1))/Nx/2.0)-1
        Nc_start = (Nx*(t_topa+1))+N_row_mid*Nx
        Nc_end=(Nx*(t_topa+1))+(N_row_mid+1)*Nx
        U_bias = Ec_old[Nc_start:Nc_end]
        Ez = (h_bar*np.pi/t_si)**2.0/(2.0*mz[0]*m_e)/q



        channel_s = junction_l
        channel_e = junction_r
        Fn_slope = (-Vd-(-Vs))/(channel_e-channel_s)

        for i_node in np.arange(0,Nx):
            if i_node < channel_s:
                Fn[i_node] = -Vs-Ez
            elif i_node >= channel_e - 1:
                Fn[i_node] = -Vd-Ez
            else:
                Fn[i_node] = -Vs-Ez+(i_node-channel_s + 1)*Fn_slope
            if fermi_flag == 1:
                Ns[i_node] = Ncc*np.log(1+np.exp((Fn[i_node]-U_bias[i_node])/(k_B*Temp/q)))
            elif fermi_flag == 0:
                Ns[i_node] = Ncc*np.exp((Fn[i_node]-U_bias[i_node])/(k_B*Temp/q))




        for iii_row in np.arange((Nx*(t_topa+1))/Nx, (Ntotal-Nx*t_bota)/Nx-1):
            for iii_col in np.arange(0, Nx):
                i_node = iii_row*Nx+iii_col
                Ne_new[i_node] = (np.sin((iii_row+1-(Nx*(t_topa+1))/Nx) * dy/t_si*np.pi)**2) * Ns[iii_col]*2.0/t_si


    #############################################################################
    elif transport_model==2: # DRIFT DIFFUSION ********************************
    #############################################################################
        mobility = np.zeros((nu_scatter-1,1))
        mobility = (np.diff(Info_scatter_old[:, 3]) + 2.0*Info_scatter_old[0:Nx-1,3])/2.0

        [U_sub, W_sub]=schred(Ec_old, Nx, Ny, Ntotal, mx, my, mz)

        E_sub=U_sub

        for i_val in np.arange(0, t_vall):
            Nccc = 2*m_e*np.sqrt(mx[i_val]*my[i_val])*k_B*Temp/(np.pi*h_bar**2)
            for i_sub in np.arange(0, max_subband):
                U_bias = U_sub[i_val, :,i_sub]
                Ns = Ne_sub_old[i_val, :,i_sub]
                if fermi_flag == 1:
                    Ns_con = Nccc*np.log(1+np.exp((-Vs-U_bias[0])/(k_B*Temp/q)))
                    Nd_con = Nccc*np.log(1+np.exp((-Vd-U_bias[Nx-1])/(k_B*Temp/q)))
                elif fermi_flag == 0:
                    Ns_con = Nccc*np.exp((-Vs-U_bias[0])/(k_B*Temp/q))
                    Nd_con = Nccc*np.exp((-Vd-U_bias[Nx-1])/(k_B*Temp/q))

                if fermi_flag == 0:
                    fac_deg = np.ones(Nx-1,1)
                elif fermi_flag == 1:
                    arg = anti_dummy(Ns/Nccc,0,fermi_flag)
                    Fn = U_bias+(k_B*Temp/q)*arg
                    zeta_channel = (Fn-U_bias)/(k_B*Temp/q)
                    fac_deg_tem = np.log(1+np.exp(zeta_channel))*(1+np.exp(-zeta_channel))
                    fac_deg_tem[np.where(arg < -5)]=1 # To ensure could scalability when subbands are empty
                    fac_deg=(np.diff(fac_deg_tem)+2*fac_deg_tem[0:Nx-1])/2
                V_channel=-U_bias
                E_channel=-np.diff(V_channel)/dx

                #If you want constant mobility uncomment the following line
                #mu_channel=mu_low./(1+(mu_low/Vel_sat*abs(E_channel)).^beta).^(1/beta)

                #If you want doping dependent mobility uncomment the following line
                mu_channel=mobility/(1+(mobility/Vel_sat*abs(E_channel))**beta)**(1/beta)
                Coe_1 = np.zeros(Nx-1)
                Coe_2 = np.zeros(Nx-1)

                for j in np.arange(0,Nx-1):
                    if abs(E_channel[j]) <= abs(e_field_limit*fac_deg[j]):
                        Coe_1[j] = -mu_channel[j]*(k_B*Temp/q)*fac_deg[j]/dx
                        Coe_2[j] = -Coe_1[j]
                    else:
                        Coe_1[j] = mu_channel[j]*E_channel[j]*1/(1-np.exp(E_channel[j]*dx/((k_B*Temp/q)*fac_deg[j])))
                        Coe_2[j] = mu_channel[j]*E_channel[j]*1/(1-np.exp(-E_channel[j]*dx/((k_B*Temp/q)*fac_deg[j])))

                AAA = np.diag(-Coe_1[1:Nx-1]+Coe_2[0:Nx-2])-np.diag(Coe_2[1:Nx-2],1)+np.diag(Coe_1[1:Nx-2],-1)
                CCC = np.zeros((Nx-2,1))
                CCC[0] = -Coe_1[0]*Ns_con/Ncc
                CCC[Nx-3] = Coe_2[Nx-2]*Nd_con/Ncc

                BBB = spsolve(sparse.csr_matrix(AAA),sparse.csr_matrix(CCC))
                #Ne_sub[i_val, :, i_sub] = np.reshape([Ns_con,BBB*Ncc,Nd_con],(Nx,1))
                Ne_sub[i_val, 0, i_sub] = Ns_con
                Ne_sub[i_val, 1:Nx-1, i_sub] = BBB*Ncc
                Ne_sub[i_val, Nx-1, i_sub] = Nd_con

                for i_node in np.arange(0,Nx):
                    N_body[:,i_node] = Ne_sub[i_val, i_node, i_sub]*W_sub[i_sub, i_val, :, i_node]/dy

                N_body_sum = N_body_sum+N_body


        Info_scatter_new[:,2] = Fn

        if ox_pnt_flag == 0:
            #Ne_new = [Ne_old(1:(Nx*(t_topa+1)));...
            #reshape((N_body_sum(2:Np_v-1,:))',...
            #Nx*(Np_v-2),1);...
            #Ne_old((Ntotal-Nx*(t_bota+1)+1):Ntotal)]
            Ne_new[0:Nx*(t_topa+1)] = Ne_old[0:Nx*(t_topa+1)]
            Ne_new[Nx*(t_topa+1):Nx*(t_topa+t_sia)] = np.reshape((N_body_sum[1:Np_v-1,:]),(1,Nx*(Np_v-2))).transpose()
            Ne_new[Nx*(t_topa+t_sia):Ntotal] = Ne_old[(Ntotal-Nx*(t_bota+1)):Ntotal]
            Ne_new = np.reshape(Ne_new, (Ntotal, 1))
        elif ox_pnt_flag == 1:
            Ne_new = np.reshape(N_body_sum.transpose(),(1, Ntotal)).transpose()


   ################################################################################
    elif transport_model == 3: #BALLISTIC TRANSPORT USING SEMICLASSICAL APPROACH#####
    ################################################################################

        [U_sub, W_sub] = schred(Ec_old, Nx, Ny, Ntotal, mx, my, mz)
        E_sub = U_sub

        for i_val in np.arange(0, t_vall):
            Ne_2d_1 = 2*(np.sqrt(mx[i_val]*my[i_val])*m_e*k_B*Temp)/(2.0*np.pi*h_bar**2)
            Ne_2d_2 = Ne_2d_1*2.0/np.pi**0.5

            for i_sub in np.arange(0,max_subband):
                U_bias = U_sub[i_val, :, i_sub]
                [Ec_peak,i_peak] = np.max(U_bias), np.argmax(U_bias)
                for i_node in np.arange (0,Nx):
                    MEc_peak = (Ec_peak-U_bias[i_node])/(k_B*Temp/q)
                    zeta_s = (-Vs-U_bias[i_node])/(k_B*Temp/q)
                    zeta_d = (-Vd-U_bias[i_node])/(k_B*Temp/q)
                    if zeta_s == zeta_d+10000:
                        Ne_sub[i_val, i_node, i_sub]=2*Ne_2d_1 * fermi(zeta_s, fermi_flag, 0)
                    else:
                        if i_node <= i_peak:
                            Ne_sub[i_val, i_node, i_sub] = Ne_2d_2*integral(zeta_s, zeta_d, MEc_peak, fermi_flag)
                        elif i_node>i_peak:
                            Ne_sub[i_val, i_node, i_sub] = Ne_2d_2*integral(zeta_d,zeta_s,MEc_peak,fermi_flag)
                    N_body[:,i_node] = Ne_sub[i_val, i_node, i_sub]*W_sub [i_sub, i_val, :,i_node]/dy
                N_body_sum = N_body_sum+N_body

        if ox_pnt_flag == 0:
            #Ne_new = np.array([Ne_old[0:(Nx*(t_topa+1))], [np.reshape((N_body_sum[1:Np_v-1,:]).transpose(),(1,Nx*(Np_v-2))).transpose()], [Ne_old[(Ntotal-Nx*(t_bota+1)):Ntotal]]])
            #Ne_old = Ne_old.toarray()
            Ne_new[0:Nx*(t_topa+1)] = Ne_old[0:Nx*(t_topa+1)]
            Ne_new[Nx*(t_topa+1):Nx*(t_topa+t_sia)] = np.reshape((N_body_sum[1:Np_v-1,:]),(1,Nx*(Np_v-2))).transpose()
            Ne_new[Nx*(t_topa+t_sia):Ntotal] = Ne_old[(Ntotal-Nx*(t_bota+1)):Ntotal]
            Ne_new = np.reshape(Ne_new,(Ntotal,1))
        elif ox_pnt_flag == 1:
            Ne_new = np.reshape(N_body_sum.transpose(), (1, Ntotal)).transpose()
    ##########################################################################
    elif transport_model == 4: #BALLISTIC TRANSPORT MODEL USING GREEN FUNCTION APPROACH
    ##########################################################################

        [U_sub , W_sub] = schred(Ec_old, Nx, Ny, Ntotal, mx, my, mz)

        E_sub = U_sub

        #E , E_number are defined outside the loop to remain constant
        # When simulating several valleys and subband so the size of the DOS matrix
        # Remain constant...This maybe more PCU intensive...

        [Ec_peak, i_peak] = np.max(U_sub[t_vall-1, :, max_subband-1]), np.argmax(U_sub[t_vall-1, :, max_subband-1])
        Ec_peak = np.max((-Vs, Ec_peak))
        E_number = round((Ec_peak+Ef_tail_up-U_sub[0, Nx-1, 0]+Ef_tail_low)/E_step)+2
        print 'E_number = ', E_number
        E = np.linspace((U_sub[0, Nx-1, 0]-Ef_tail_low), (Ec_peak+Ef_tail_up), E_number)
        delta_E = ((Ec_peak+Ef_tail_up)-(U_sub[0, Nx-1, 0]-Ef_tail_low))/(E_number-1)

        for i_val in np.arange(0, t_vall):
            Ne_2d = 2*np.sqrt(mx[i_val]*m_e*(k_B*Temp/q)*q/(2*np.pi**3))/(h_bar*dx)
            tt = (h_bar**2)/(2*my[i_val]*m_e*(dx**2)*q)
            tt = float(tt)
            A = tt*((2*np.eye(Nx))-(np.diag(np.ones(Nx-1), 1))-(np.diag(np.ones(Nx-1), -1)))

            for i_sub in np.arange(0,max_subband):
                U_bias = U_sub[i_val, :, i_sub]
                #[Ec_peak,i_peak]=max(U_bias)
                #Ec_peak=max(-Vs,max(U_bias))
                #E_number=round((Ec_peak+Ef_tail_up-U_bias(Nx)+Ef_tail_low)/E_step)+2
                #E=linspace((U_bias(Nx)-Ef_tail_low),(Ec_peak+Ef_tail_up),E_number)
                #delta_E=((Ec_peak+Ef_tail_up)-(U_bias(Nx)-Ef_tail_low))/(E_number-1)
                N_den = np.zeros((Nx,1))
                B_s = np.zeros((Nx,1))
                B_d = np.zeros((Nx,1))
                B_s[0] = 1
                B_d[Nx-1] = 1
                spB_s = sparse.csr_matrix(B_s)
                spB_d = sparse.csr_matrix(B_d)

                #for k=1:E_number,
                    #ee=E(k);ep=ee+eta
                    #ck=1-((ep-U_bias(1))/(2*tt));con_s=-tt*exp(i*acos(ck))
                    #ck=1-((ep-U_bias(Nx))/(2*tt));con_d=-tt*exp(i*acos(ck))
                    #U_eff=U_bias
                    #U_eff(1)=U_bias(1)+con_s
                    #U_eff(Nx)=U_bias(Nx)+con_d
                    #G_inv=sparse((ep*eye(Nx))-A-diag(U_eff))
                    #G_s=G_inv\spB_s
                    #G_d=G_inv\spB_d
                    #f_1=fermi(((-Vs-ee)/(k_B*Temp/q)),fermi_flag,-1/2)
                    #f_2=fermi(((-Vd-ee)/(k_B*Temp/q)),fermi_flag,-1/2)
                    #N_den=N_den-abs(G_s).^2*imag(con_s)*2*f_1...
                    #      -abs(G_d).^2*imag(con_d)*2*f_2
                    #N_dos_one(:,k)=-abs(G_s).^2*imag(con_s)-abs(G_d).^2*imag(con_d)
                #end

                [N_den1, Nquad] = myquad(func_energy, E[0], E[E_number-1], 1e-6, [], tt, U_bias, A, spB_s, spB_d)
                N_den = N_den1/(E[1]-E[0])

                Ne_sub[i_val, :, i_sub] = (N_den)*Ne_2d*delta_E
                for i_node in np.arange(0,Nx):
                    N_body[:,i_node] = Ne_sub[i_val, i_node, i_sub]*W_sub[i_sub, i_val, :, i_node]/dy

                N_body_sum = N_body_sum+N_body

        if ox_pnt_flag == 0:
            #Ne_new=[Ne_old[1:(Nx*(t_topa+1)));...
            #reshape((N_body_sum(2:Np_v-1,:))',...
            #Nx*(Np_v-2),1);...
            #Ne_old((Ntotal-Nx*(t_bota+1)+1):Ntotal)]
            Ne_new[0:Nx*(t_topa+1)] = Ne_old[0:Nx*(t_topa+1)]
            Ne_new[Nx*(t_topa+1):Nx*(t_topa+t_sia)] = np.reshape((N_body_sum[1:Np_v-1,:]),(1,Nx*(Np_v-2))).transpose()
            Ne_new[Nx*(t_topa+t_sia):Ntotal] = Ne_old[(Ntotal-Nx*(t_bota+1)):Ntotal]
            Ne_new = np.reshape(Ne_new, (Ntotal, 1))
        elif ox_pnt_flag == 1:
           Ne_new = np.reshape(N_body_sum.transpose(), (1, Ntotal)).transpose()

        globvars.E = E

    ################################################################################
    elif transport_model == 5: #SCATTERING MODEL USING GREEN"S FUNCTION METHOD*******
    ################################################################################

        [U_sub, W_sub] = schred(Ec_old, Nx, Ny, Ntotal, mx, my, mz)
        E_sub = U_sub

        U_scatter = np.zeros((Nx, 1))
        mu_scatter_old = np.zeros((nu_scatter+2, 1))
        mu_scatter_new = np.zeros((nu_scatter+2, 1))
        delta_mu = np.zeros((nu_scatter+2, 1))
        i_scatter = np.zeros((nu_scatter+2, 1))
        BB_dummy = np.eye(Nx)
        Iin = np.zeros((nu_scatter,1))

        #comment the two following line for constant mobility
        #zeta_self=zeros(nu_scatter,1);
        #decay_fac=1;%dimensionless

        for i_s in np.arange(0, nu_scatter):
            i_scatter[i_s] = (Info_scatter_old[i_s, 1])

            # CHANGED BY RAMESH TO TAKE OLD GUESS
            mu_scatter_old[i_s] = (Info_scatter_new[i_s, 2])

        mu_scatter_old[nu_scatter] = -Vs
        mu_scatter_old[nu_scatter+1] = -Vd
        i_scatter[nu_scatter] = 0
        i_scatter[nu_scatter+1] = Nx-1
        BB = sparse.lil_matrix((int(Nx), int(nu_scatter+2)), dtype=float)
        i_scatter = i_scatter.astype(int)
        i_scatter = np.squeeze(i_scatter)
        #mu_scatter_old = mu_scatter_old.astype(int)
        mu_scatter_old = np.squeeze(mu_scatter_old)

        BB[:, 0:int(nu_scatter)] = BB_dummy[:,i_scatter[0:nu_scatter]]
        BB[:, int(nu_scatter):int(nu_scatter)+1] = BB_dummy[:, 0:1]
        for ind in np.arange(0,Nx):
            if BB_dummy[ind, Nx-1]:
                BB[ind, int(nu_scatter) + 1] = BB_dummy[ind, Nx-1]
        BB = BB.tocsr()

        Ec_peak = np.max((-Vs, np.max(U_sub[0, :, 0])))
        Ec_bottom = U_sub[0, Nx-1, 0]
        E_number = round((Ec_peak+Ef_tail_up-Ec_bottom+Ef_tail_low)/E_step)+2
        E = np.linspace((Ec_bottom-Ef_tail_low), (Ec_peak+Ef_tail_up), E_number)
        delta_E = ((Ec_peak+Ef_tail_up)-(Ec_bottom-Ef_tail_low))/(E_number-1)

        GG = sparse.csr_matrix((Nx,nu_scatter+2))
        T_E = np.zeros((nu_scatter+2, E_number, nu_scatter+2))

        # comment the next line and uncomment the previous line for constant mobility

        E_self = np.ones((E_number,nu_scatter+2)) #correction on March 14th to include mobility
        E_self = E_self + 0j
        U_tem = np.zeros((nu_scatter,1))

        # end of initialization

        # ADDED BY RAMESH JAN 13th 03.
        # We need to compute a mean free path to set
        # zeta_self based on the Fermi-level of the
        # probe, subband energy and grid spacing.
        sum1 = np.zeros((nu_scatter, 1))
        sum2 = np.zeros((nu_scatter, 1))
        lambda1 = np.zeros((nu_scatter, 1))
        for i_val in np.arange(0, t_vall):
            for i_sub in np.arange(0, max_subband):
                for i_s in np.arange(0, nu_scatter):
                    factor = np.exp((Info_scatter_new[i_s,2]-U_sub[i_val, i_scatter[i_s], i_sub])/(k_B*Temp/q))
                    sum1[i_s] = sum1[i_s]+my[i_val]*m_e*np.log(1+factor)
                        # Prevent divide by zero
                    factor1 = np.log(1+factor)
                    if(factor1<=1e-20):
                        factor1 = factor
                    factor2 = factor/(factor1*(1+factor))*fermi(((Info_scatter_new[i_s,2]-U_sub[i_val, i_scatter[i_s],i_sub])/(k_B*Temp/q)),fermi_flag,1.0/2.0)
                    sum2[i_s] = sum2[i_s]+my[i_val]*m_e*np.sqrt(2*k_B*Temp/(np.pi*mx[i_val]*m_e))*factor2

        for i_s in np.arange(0, nu_scatter):
            lambda1[i_s]=2*k_B*Temp/q*Info_scatter_old[i_s,3]*sum1[i_s]/sum2[i_s]

        lambda1 = np.squeeze(lambda1)

        ttot=0
        scattot=0
        ktot =0
        stot = 0
        k2tot =0
        soltot = 0
        asstot = 0
        N_dos = np.zeros((Nx, E_number))

        for i_val in np.arange(0,t_vall):
            t1 = time.time()
            Ie_2d = 2*q**2/(np.pi**2*h_bar**2)*np.sqrt(my[i_val]*m_e*(k_B*Temp/q)*q*np.pi/2.0)
            tt = (h_bar**2)/(2*mx[i_val]*m_e*(dx**2)*q)
            tt = float(tt)
            A = tt*((2*np.eye(Nx))-(np.diag(np.ones(Nx-1),1))-(np.diag(np.ones(Nx-1),-1)))
            t2 = time.time()
            ttot += (t2 - t1)
            for i_sub in np.arange(0, max_subband):
                s1 = time.time()
                U_bias = U_sub[i_val, :, i_sub]
                U_tem = U_bias[i_scatter[0:nu_scatter]]
                s2 = time.time()
                stot += (s2-s1)
                for k in np.arange(0, E_number):
                    k1=time.time()
                    ee = E[k]
                    ep = ee+eta
                    ck = 1-((ep-U_bias[0])/(2*tt))
                    con_s = -tt*np.exp(1j*np.arccos(ck))
                    ck = 1-((ep-U_bias[Nx-1])/(2*tt))
                    con_d = -tt*np.exp(1j*np.arccos(ck))
                    k2 = time.time()
                    ktot += (k2 -k1)
                # ADDED BY RAMESH JAN 13 03.
                    scat1 = time.time()
                    for i_s in np.arange(0, nu_scatter):
                        elutot = ee-(A[i_scatter[i_s], i_scatter[i_s]]+U_bias[i_scatter[i_s]])
                        g1 = (elutot+cmath.sqrt(elutot**2-4.0*tt**2))/(2*tt**2)
                        g2 = (elutot-cmath.sqrt(elutot**2-4.0*tt**2))/(2*tt**2)
                        g3 = (np.imag(g2)<=0)*g2+(np.imag(g1<=0))*g1
                        E_self[k,i_s] = 1j*np.imag(g3*2*dx*tt**2/lambda1[i_s])

                    scat2 = time.time()
                    scattot += (scat2-scat1)
                    k3 = time.time()
                    E_self[k, nu_scatter] = con_s
                    E_self[k, nu_scatter+1] = con_d
                    U_scatter = U_scatter + 0j
                    U_scatter = np.squeeze(U_scatter)
                    U_scatter[i_scatter] = (E_self[k, :])

                    U_eff = U_bias+U_scatter
                    #print U_eff
                    G_inv = sparse.csr_matrix((ep*np.eye(Nx))-A-np.diag(U_eff))
                    GG = np.zeros((Nx, nu_scatter+2))
                    GG = GG + 0j
                    k4 = time.time()
                    k2tot += (k4-k3)
                    sol1 = time.time()
                    gin = np.linalg.inv(G_inv.todense())
                    gin = sparse.csr_matrix(gin)
                    GG = gin*BB
                    GG = GG.toarray()
                    sol2 = time.time()
                    soltot += (sol2-sol1)
                    #print sparse.csr_matrix(np.diag(np.imag(E_self[k,:]).conj().T))
                    #print np.shape(GG[i_scatter,:])
                    #print delta_E
                    #print abs(GG[i_scatter,:])**2
                    #print sparse.csr_matrix(4*int(Ie_2d)*delta_E*np.dot(np.diag(np.imag(E_self[k,:]).conj().T),abs(GG[i_scatter,:])**2))
                    T_E[:,k,:] = 4*int(Ie_2d)*delta_E*np.dot((np.dot(np.diag(np.imag(E_self[k,:]).conj().T),abs(GG[i_scatter,:])**2)),np.diag(np.imag(E_self[k,:]).conj().T))+np.squeeze(T_E[:, k, :])
                    #print sparse.csr_matrix(T_E[:,k,:])
                    #print E_self[k,:]
                    #print GG[i_scatter,:]
                    #print k
                    sol3 = time.time()
                    asstot += (sol3-sol2)
        #print ttot, stot, ktot, scattot, k2tot, soltot, asstot

        #calculate current, and search for mu_scatter
        #based on current continuity

        [Is, Id, Iin, mu_scatter_new] = current_mat(mu_scatter_old, T_E, E)
        #print 'sad'
        #print Is
        #print Id

        #end of mu_scatter search

        #Calculate charge density
        for i_val in np.arange(0, t_vall):
            Ne_2d = 2.0*np.sqrt(my[i_val]*m_e*(k_B*Temp/q)*q/(2.0*np.pi**3))/(h_bar*dx)
            tt = (h_bar**2)/(2.0*mx[i_val]*m_e*(dx**2)*q)
            tt = float(tt)
            A = tt*((2.0*np.eye(Nx))-(np.diag(np.ones(Nx-1),1))-(np.diag(np.ones(Nx-1),-1)))

            for i_sub in np.arange(0,max_subband):
                U_bias = U_sub[i_val, :, i_sub]
                U_tem = U_bias[(i_scatter[0:nu_scatter])]

                for k in np.arange(0, E_number):

                    # All this needs to be recomputed for multiple bands
                    #--------------------------------------------------
                    ee=E[k]
                    ep=ee+eta
                    ck=1-((ep-U_bias[0])/(2*tt))
                    con_s=-tt*np.exp(1j*np.arccos(ck))
                    ck=1-((ep-U_bias[Nx-1])/(2*tt))
                    con_d=-tt*np.exp(1j*np.arccos(ck))

                 # ADDED BY RAMESH JAN 13 03.
                    for i_s in np.arange(0, nu_scatter):
                        elutot = ee-(A[i_scatter[i_s], i_scatter[i_s]]+U_bias[i_scatter[i_s]])
                        g1 = (elutot+cmath.sqrt(elutot**2-4*tt**2))/(2*tt**2)
                        g2 = (elutot-cmath.sqrt(elutot**2-4*tt**2))/(2*tt**2)
                        g3 = (np.imag(g2)<=0)*g2+(np.imag(g1<=0))*g1
                        E_self[k,i_s] = 1j*np.imag(g3*2*dx*tt**2/lambda1[i_s])

                    E_self = E_self + 0j
                    E_self[k, nu_scatter] = con_s
                    E_self[k, nu_scatter+1] = con_d
                    U_scatter[i_scatter] = (E_self[k,:]).transpose()

                    U_eff = U_bias+U_scatter
                    G_inv = sparse.csr_matrix((ep*np.eye(Nx))-A-np.diag(U_eff))
                    gin = np.linalg.inv(G_inv.todense())
                    gin = sparse.csr_matrix(gin)
                    GG = gin*BB
                    GG = GG.toarray()
                    Ne_sub[i_val, :, i_sub]=Ne_sub[i_val, :, i_sub]-np.dot(2*Ne_2d*delta_E*abs(GG)**2,(np.imag(E_self[k,:]).conj().T)*fermi(((mu_scatter_new-ee)/(k_B*Temp/q)),fermi_flag,-1.0/2.0))

                    GGG = np.zeros((Nx,Nx))
                    GGG = GGG +0j
                    GGG[:, 0] = GG[:, Nx-2]
                    GGG[:, Nx-1] = GG[:, Nx-1]
                    GGG[0:Nx, 1:Nx-1] = GG[0:Nx,0:Nx-2]
                    N_dos[:,k] = -2*np.imag(np.diag(GGG))


                for i_node in np.arange(0,Nx):
                    N_body[:, i_node] = Ne_sub[i_val, i_node, i_sub]*W_sub[i_sub, i_val, :,i_node]/dy

                N_body_sum = N_body_sum+N_body

        Info_scatter_new[:, 0] = np.squeeze(Iin)
        Info_scatter_new[0, 0] = Is
        Info_scatter_new[:, 2] = mu_scatter_new[0:nu_scatter]

        if ox_pnt_flag==0:
            #Ne_new=[Ne_old(1:(Nx*(t_topa+1)));...
            #reshape((N_body_sum(2:Np_v-1,:))',...
            #Nx*(Np_v-2),1);...
            #Ne_old((Ntotal-Nx*(t_bota+1)+1):Ntotal)]
            Ne_new[0:Nx*(t_topa+1)] = Ne_old[0:Nx*(t_topa+1)]
            Ne_new[Nx*(t_topa+1):Nx*(t_topa+t_sia)] = np.reshape((N_body_sum[1:Np_v-1,:]),(1,Nx*(Np_v-2))).transpose()
            Ne_new[Nx*(t_topa+t_sia):Ntotal] = Ne_old[(Ntotal-Nx*(t_bota+1)):Ntotal]
            Ne_new = np.reshape(Ne_new, (Ntotal, 1))
        elif ox_pnt_flag==1:
            Ne_new=np.reshape(N_body_sum.transpose(), (1, Ntotal)).transpose()

        globvars.E = E
        globvars.Info_scatter_new = Info_scatter_new
        globvars.Info_scatter_old = Info_scatter_old
        globvars.Is = Is
        globvars.Id = Id
        globvars.N_dos = N_dos


    ################################################################################
    ######################START OF VARIABLE CHANGE PART#############################
    ################################################################################

    if ox_pnt_flag == 0:  # No electron penetration into oxide
        for iii_row in np.arange((Nx*(t_topa+1))/Nx, (Ntotal-Nx*t_bota)/Nx-1):
            for iii_col in np.arange(0, Nx):
                i_node = iii_row*Nx+iii_col
                Fn_new[i_node] = anti_dummy(Ne_new[i_node]/Nc, dummy_flag, fermi_flag)*(k_B*Temp/q)+Ec_old[i_node]
    elif ox_pnt_flag == 1: # assuming electron penetration into oxide
        Fn_new = anti_dummy((Ne_new+div_avd)/Nc, dummy_flag, fermi_flag)*(k_B*Temp/q)+Ec_old

    #############################################################################
    #######################END OF VARIABLE CHANGE PART###########################
    #############################################################################

    #############################################################################
    ########################END OF FUNCTION CHARGE.M#############################
    #############################################################################

    return [Fn_new, Ne_new, Ne_sub, E_sub]