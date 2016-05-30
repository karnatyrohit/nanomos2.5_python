#################################################################################
#CALCULATE THE CURRENT DENSITY [A/m] GIVEN A 1D CHARGE SHEET POTENTIAL PROFILE###
#(Zhibin Ren 7-18-00)############################################################
#################################################################################

from readinput import *
from fermi import fermi

def current(Ne, Ec, NE_sub, E_sub, Nx, Ny, Ntotal, mx, my ,mz ):
    Temp = Te
    transport_model = transportmodel.value
    Vd = Vdc.value
    fermi_flag = fermiflag1.value

    #### MODEL 2 related parameters

    eta = 1e-6j
    Ef_tail_low = 0.01
    Ef_tail_up = 0.3
    E_step = criterion_outer/2  # 0.5 times criterion_outer

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
    Ie_sub = np.zeros((Nx,max_subband,t_vall))
    Te_sub = np.zeros((Nx,max_subband,t_vall))
    Fn_sub = np.zeros((Nx,max_subband,t_vall))
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

    ###########################################################################
    elif transport_model==3:  #BALLISTIC TRANSPORT USING SEMICLASSICAL APPROACH
    ###########################################################################
        U_bias = np.zeros((Nx,1))

        for i_val in np.arange(0, t_vall):
            Ie_2d = 2*q/h_bar**2*np.sqrt(mx[i_val]*m_e/2)*((k_B*Temp/q)*q/np.pi)**(3/2)

            for i_sub in np.arange(0, max_subband):
                U_bias = E_sub[i_val, :,i_sub]
                Ec_peak = max(U_bias)

                Ie_tem = 0
                Ie_tem = Ie_2d*(fermi(((-Vs-Ec_peak)/(k_B*Temp/q)),fermi_flag,1/2)-fermi(((-Vd-Ec_peak)/(k_B*Temp/q)),fermi_flag,1/2))
                Ie_sub[i_val, :, i_sub] = Ie_tem*np.ones(Nx)
                Ie = Ie+Ie_tem*np.ones(Nx)

    return [Ie, Ie_sub, Te_sub, Fn_sub]

