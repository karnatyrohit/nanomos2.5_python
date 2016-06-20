#################################################################################
#CALCULATE THE CURRENT DENSITY [A/m] GIVEN A 1D CHARGE SHEET POTENTIAL PROFILE###
#(Zhibin Ren 7-18-00)############################################################
#################################################################################

from readinput import *
from fermi import fermi
from globvars import globvars
from scipy import  sparse
from scipy.sparse.linalg import spsolve
from anti_dummy import anti_dummy

def current(Ne, Ec, Ne_sub, E_sub, Nx, Ny, Ntotal, mx, my ,mz ):
    Temp = Te
    transport_model = transportmodel.value
    Vd = Vdc.value
    fermi_flag = fermiflag1.value
    E = globvars.E
    nu_scatter = globvars.nu_scatter
    Info_scatter_old = globvars.Info_scatter_old

    #### MODEL 2 related parameters

    eta = 1e-6j
    Ef_tail_low = 0.01
    Ef_tail_up = 0.3
    E_step = criterion_outer/2.0  # 0.5 times criterion_outer

    ####MODEL 4 related parameters
    e_field_limit = 1.0  # in V/m
    ####MODEL 5 related parameters

    ###########################################################
    #	INPUT AND OUTPUT VARIABLES
    ###########################################################
    #Info_scatter_old contains all information of scatters
    #2-dim matrix, dim_2:current,index, and Fermi energy, low field mobility
    #Ne_old is the old electron density
    #1 column of Ntotal elements
    #Ec_old is the old potential energy profile in eV
    #1 column of Ntotal elements
    #Te_old is the old electron temperature along the channel
    #1 column of Nx elements, for charge-sheet-model
    #Ne_sub contains 2D electron density on each subband
    #E_sub contains the 1D subband energy for all subband
    #considered
    #Te_sub_old contains the 1D electron temperature for all subband
    #considered, for bulk model

    #Ne is the 3D electron density profile in the entire device
    #1 column matrix of Ntotal elements
    #Ec is the conduction band edge profile in the entire device
    #1 column matrix of Ntotal elements
    #Te_old is 1 column matrix of Nx elements for electron temperature
    #, for charge-sheet-model

    #Ie is 1 column matrix of Nx elements for current density
    #A/m
    #Ie_sub contains currents carried in all subbands.
    ###########################################################

    ###########################################################
    #	TEMPORARY VARIABLES
    ###########################################################
    #Ec_channel is the given potential energy profile
    #along the channel
    #Fn_channel is the given quasi Fermi level along the channel
    #Ne_channel is the given electron 2D density along the channel
    #max_subband is the number of subbands considered
    #t_vall is the number of valleyes considered

    #########################INITIALIZATION#######################################
    Ie = np.zeros((Nx,1))
    Ie_sub = np.zeros((t_vall, Nx,max_subband))
    Te_sub = np.zeros((t_vall, Nx,max_subband))
    Fn_sub = np.zeros((t_vall, Nx,max_subband))
    #E_sub = np.zeros((Nx,max_subband,t_vall))
    #Ne_sub = np.zeros((Nx,max_subband,t_vall))

    #Info_scatter_new = np.zeros((nu_scatter,4))
    #Info_scatter_new = Info_scatter_old
    Ec_channel = np.zeros((Nx,1))
    Fn_channel = np.zeros((Nx,1))
    Ne_channel = np.zeros((Nx,1))
    #Fn_new = np.zeros((Ntotal,1))
    #Ne_new = np.zeros((Ntotal,1))

    #temporarily used variables
    if ox_pnt_flag == 0:
       Np_v = round(t_si/dy)+1
    elif ox_pnt_flag == 1:
       Np_v = Ny

    Ns = np.zeros((Nx, 1))
    N_body = np.zeros((Np_v, Nx))
    N_body_sum = np.zeros((Np_v, Nx))
    U_sub = np.zeros((Nx, max_subband, t_vall))
    W_sub = np.zeros((Np_v, Nx, max_subband, t_vall))
    U_bias = np.zeros((Nx, 1))

    ##############################################################################
    #THE START OF TRANSPORT EQUATION PART#########################################
    ##############################################################################

    ##############################################################################
    if transport_model==1:  #TRANSPORT MODEL 1######################################
        print 'entered tm1'
    ##############################################################################

    #*****************************************************************************
    elif transport_model==2: #TRANSPORT MODEL 2 DD*******************************
    #*****************************************************************************
        U_bias = np.zeros((Nx,1))
        Ns = np.zeros((Nx,1))
        Ie_tem = np.zeros((Nx,1))
        mobility = np.zeros((nu_scatter-1,1))
        mobility = (np.diff(Info_scatter_old[:,3])+2*Info_scatter_old[0:Nx-1,3])/2.0

        for i_val in np.arange(0, t_vall):
            Nccc = 2*m_e*np.sqrt(mx[i_val]*my[i_val])*k_B*Temp/(np.pi*h_bar**2)

            for i_sub in np.arange(max_subband):
                U_bias = E_sub[i_val, :, i_sub]
                Ns = Ne_sub[i_val, :, i_sub]

                if fermi_flag == 0:
                    fac_deg= np.ones(Nx-1)
                elif fermi_flag == 1:
                    Fn_channel = U_bias+(k_B*Temp/q)*anti_dummy(Ns/Nccc,0,fermi_flag)
                    zeta_channel = (Fn_channel-U_bias)/(k_B*Temp/q)
                    fac_deg_tem = np.log(1+np.exp(zeta_channel))*(1+np.exp(-zeta_channel))
                    fac_deg=(np.diff(fac_deg_tem)+2*fac_deg_tem[0:Nx-1])/2
                V_channel = -U_bias
                E_channel = -np.diff(V_channel)/dx
                mu_channel = mobility/(1+(mobility/Vel_sat*abs(E_channel))**beta)**(1/beta)

                Coe_1 = np.zeros(Nx-1)
                Coe_2 = np.zeros(Nx-1)

                for j in np.arange(0,Nx-1):
                    if abs(E_channel[j])<=abs(e_field_limit*fac_deg[j]):
                        Coe_1[j] = -mu_channel[j]*(k_B*Temp/q)*fac_deg[j]/dx
                        Coe_2[j] = -Coe_1[j]
                    else:
                        Coe_1[j] = mu_channel[j]*E_channel[j]*1/(1-np.exp(E_channel[j]*dx/((k_B*Temp/q)*fac_deg[j])))
                        Coe_2[j] = mu_channel[j]*E_channel[j]*1/(1-np.exp(-E_channel[j]*dx/((k_B*Temp/q)*fac_deg[j])))
                    Ie_tem[j]=-q*(Ns[j]*Coe_1[j]+Ns[j+1]*Coe_2[j])

                Ie_tem[Nx-1] = Ie_tem[Nx-2]
                Ie_sub[i_val, :,i_sub] = np.squeeze(Ie_tem)
                Fn_sub[i_val, :,i_sub] = Fn_channel
                Ie = Ie+Ie_tem

    ###########################################################################
    elif transport_model==3:  #BALLISTIC TRANSPORT USING SEMICLASSICAL APPROACH
    ###########################################################################
        U_bias = np.zeros((Nx,1))

        for i_val in np.arange(0, t_vall):
            Ie_2d = 2.0*q/h_bar**2*np.sqrt(mx[i_val]*m_e/2.0)*((k_B*Temp/q)*q/np.pi)**(3.0/2.0)

            for i_sub in np.arange(0, max_subband):
                U_bias = E_sub[i_val, :,i_sub]
                Ec_peak = max(U_bias)

                Ie_tem = 0
                Ie_tem = Ie_2d*(fermi(((-Vs-Ec_peak)/(k_B*Temp/q)),fermi_flag,1.0/2.0)-fermi(((-Vd-Ec_peak)/(k_B*Temp/q)),fermi_flag,1.0/2.0))
                Ie_sub[i_val, :, i_sub] = Ie_tem*np.ones(Nx)
                Ie = Ie+Ie_tem*np.ones(Nx)


    ##################################################################################
    elif transport_model==4: #BALLISTIC TRANSPORT MODEL USING GREEN FUNCTION APPROACH#
    ##################################################################################

        U_bias = np.zeros((Nx,1))

        Ec_peak = np.max((-Vs, np.max(E_sub[t_vall-1, :, max_subband-1])))
        E_number = round((Ec_peak+Ef_tail_up - E_sub[0, Nx-1, 0] + Ef_tail_low)/E_step)+2
        E = np.linspace((E_sub[0, Nx-1, 0]-Ef_tail_low),(Ec_peak+Ef_tail_up),E_number)
        delta_E = ((Ec_peak+Ef_tail_up)-(E_sub[0, Nx-1, 0]-Ef_tail_low))/(E_number-1)

        N_dos_one = np.zeros((Nx,E_number))
        N_dos = np.zeros((Nx,E_number))

        #used to plot the transmission coefficient with several valleys
        Trans_sub = np.zeros((E_number,1))
        Trans = np.zeros((E_number,1))

        for i_val in np.arange(0,t_vall):
            Ne_2d = 2*np.sqrt(mx[i_val]*m_e*(k_B*Temp/q)*q/(2*np.pi**3))/(h_bar*dx)
            Ie_2d = 2*q**2/(np.pi**2*h_bar**2)*np.sqrt(mx[i_val]*m_e*(k_B*Temp/q)*q*np.pi/2)
            tt = (h_bar**2)/(2*my[i_val]*m_e*(dx**2)*q)
            tt = float(tt)
            A = tt*((2*np.eye(Nx))-(np.diag(np.ones(Nx-1), 1))-(np.diag(np.ones(Nx-1), -1)))

            for i_sub in np.arange(0,max_subband):
                U_bias = E_sub[i_val, :, i_sub]
                #Ec_peak = max(-Vs,max(U_bias))
                #E_number = round((Ec_peak+Ef_tail_up-U_bias(Nx)+Ef_tail_low)/E_step)+2
                #E = np.linspace((U_bias(Nx)-Ef_tail_low),(Ec_peak+Ef_tail_up),E_number)
                #delta_E = ((Ec_peak+Ef_tail_up)-(U_bias(Nx)-Ef_tail_low))/(E_number-1)

                B_d = np.zeros((Nx, 1))
                B_d[Nx-1] = 1
                spB_d = sparse.csr_matrix(B_d)
                B_s = np.zeros((Nx, 1))
                B_s[0] = 1
                spB_s = sparse.csr_matrix(B_s)

                Ie_tem = 0
                for k in np.arange(0,E_number):
                    ee = E[k]
                    ep = ee + eta
                    ck = 1-((ep-U_bias[0])/(2*tt))
                    con_s = -tt*np.exp(1j*np.arccos(ck))
                    ck = 1-((ep-U_bias[Nx-1])/(2*tt))
                    con_d = -tt*np.exp(1j*np.arccos(ck))
                    U_eff = U_bias
                    U_eff = U_eff + 0j
                    U_eff[0] = U_bias[0] + con_s
                    U_eff[Nx-1] = U_bias[Nx-1] + con_d
                    G_inv = sparse.csr_matrix((ep*np.eye(Nx))-A-np.diag(U_eff))
                    G_s = spsolve(G_inv, spB_s)
                    G_d = spsolve(G_inv, spB_d)
                    f_1 = fermi(((-Vs-ee)/(k_B*Temp/q)), fermi_flag, -1.0/2.0)
                    f_2 = fermi(((-Vd-ee)/(k_B*Temp/q)), fermi_flag, -1.0/2.0)
                    Ie_tem = Ie_tem+abs(G_d[0])**2*np.imag(con_s)*np.imag(con_d)*4*(f_1-f_2)
                    Trans[k] = abs(G_d[0])**2*np.imag(con_s)*np.imag(con_d)*4
                    N_dos_one[:, k] = -abs(G_s)**2*np.imag(con_s)-abs(G_d)**2*np.imag(con_d)

                N_dos = N_dos+Ne_2d*N_dos_one

                Ie_sub[i_val, :, i_sub] = Ie_2d*delta_E*Ie_tem * np.ones(Nx)
                Ie = Ie+Ie_2d*delta_E * Ie_tem*np.ones(Nx)
                Trans_sub = Trans_sub+Trans

        Trans = Trans_sub
        globvars.Trans = Trans
        globvars.N_dos = N_dos

    return [Ie, Ie_sub, Te_sub, Fn_sub]

