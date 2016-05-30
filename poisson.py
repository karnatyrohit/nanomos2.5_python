###########################################################################
#########2D Poisson equation solver using Newton approach (Zhibin Ren 7-18-00)
###########################################################################

from scipy import sparse
from readinput import *
from dummy import dummy
from dummy_prime import dummy_prime

def poisson(spNd, spFn, Ec_old, F_prime, div_avd, charge_fac, Eg1, Eg2, Es, Ed, Nx, Ny, Ntotal):
    Temp = Te
    transport_model = transportmodel.value
    fermi_flag = fermiflag1.value

    Lsda = round(Lsd/dx)
    Lg_topa = round(Lg_top/dx)
    Lg_bota = round(Lg_bot/dx)
    t_topa = round(t_top/dy)
    t_bota = round(t_bot/dy)
    t_sia = round(t_si/dy)

    ###########################INPUT AND OUTPUT VARIABLES######################
    # spNd is the 3d net dopant density
    #   1 column of Ntotal element
    # spFn is the 3d electron density converted quasi Fermi energy
    #   from the solutions of transport equation
    #   1 column of Ntotal elements
    # Ec_old is the input old potential energy profile
    #   1 column of Ntotal elements
    # F_prime is the derivative matrix of Ne=0 F (Ntotal by Ntotal)
    # criterion_inner is the poisson solver convergence criterion in eV
    # Ec is the new potential profile computed
    #   1 column of Ntotal elements
    # All equations involved assume ISS UNITS <- KDC No such thing!
    # Uses SI units (kg-m-s) I think

    ###############################INITIALIZATION##############################
    delta_Ec = np.zeros((Ntotal,1))
    F = np.zeros((Ntotal,1))
    MF_prime = np.zeros((Ntotal,Ntotal))
    Ec = np.zeros((Ntotal,1))
    Fn = np.zeros((Ntotal,1))
    Nd = np.zeros((Ntotal,1))
    dummy_fun = np.zeros((Ntotal,1))
    dummy_fun_prime = np.zeros((Ntotal,1))

    iter_inner = 0
    error_inner = 1
    Ec = np.real(Ec_old)
    Nd = np.real(spNd)
    Fn = np.real(spFn)

    if ox_pnt_flag == 0:
        div_avdt = 0
    elif ox_pnt_flag == 1:
        div_avdt = div_avd

    ##########################START OF INNER LOOP######################################

    while error_inner >= criterion_inner:
        iter_inner = iter_inner+1
        print '%s %i \n' % ('iter_inner = ',iter_inner)

        ####################THE START OF DUMMY VARIABLE DEFINITION#################
        dummy_fun = charge_fac * ((Nd+div_avdt)/Nc - dummy((Fn-Ec)/(k_B*Temp/q), dummy_flag,fermi_flag))
        dummy_fun_prime = charge_fac * dummy_prime((Fn-Ec)/(k_B*Temp/q), dummy_flag,fermi_flag)/(k_B*Temp/q)


        ########################THE END OF DUMMY VARIABLE DEFINITION###############

        if ox_pnt_flag == 0: # NO ELECTRON PENETRATION INTO OXIDE REGIONS

    #################################EVALUATE F#########################################

    #############################Top gate insulator region##############################
            for i in range(0,Nx*(t_topa+1)):
                if(i >= 0 and i < Lsda):
                    F[i] = Ec[i]-Ec[i+Nx]
                elif(i >= Lsda and i <= (Lsda+Lg_topa)):
                    F[i] = Ec[i]-Eg1
                elif(i >= (Lsda+Lg_topa)+1 and i<Nx):
                    F[i] = Ec[i]-Ec[i+Nx]
                elif(i >= (Nx*t_topa) and i < Nx*(t_topa+1)):
                    F[i] = -1/8*(dy/dx)*eps_top/eps_si*Ec[i-Nx-1] \
                        - (dx/dy-1/4/(dx/dy))*eps_top/eps_si*Ec[i-Nx] \
                        - 1/8*(dy/dx)*eps_top/eps_si*Ec[i-Nx+1] \
                        - 3/8*(dy/dx)*(eps_top+eps_si)/eps_si*Ec[i-1] \
                        + (dx/dy+3/4/(dx/dy))*(eps_top+eps_si)/eps_si*Ec[i] \
                        - 3/8*(dy/dx)*(eps_top+eps_si)/eps_si*Ec[i+1] \
                        - 1/8*(dy/dx)*Ec[i+Nx-1] \
                        - (dx/dy-1/4/(dx/dy))*Ec[i+Nx] \
                        - 1/8*(dy/dx)*Ec[i+Nx+1]
                else:
                    F[i] = -(dx/dy)*eps_top/eps_si*Ec[i-Nx] \
                        -(dy/dx)*eps_top/eps_si*Ec[i-1] \
                        +2*(dx/dy+dy/dx)*eps_top/eps_si*Ec[i] \
                        -(dy/dx)*eps_top/eps_si*Ec[i+1] \
                        -(dx/dy)*eps_top/eps_si*Ec[i+Nx]

            #########################Bottom gate insulator region##############################
            for i in range((Ntotal-Nx*(t_bota+1)),Ntotal):
                if(i >= (Ntotal-Nx*(t_bota+1)) and i<(Ntotal-Nx*t_bota)):
                    F[i] = -1/8*(dy/dx)*Ec[i-Nx-1] \
                        - (dx/dy-1/4/(dx/dy))*Ec[i-Nx] \
                        - 1/8*(dy/dx)*Ec[i-Nx+1] \
                        - 3/8*(dy/dx)*(eps_bot+eps_si)/eps_si*Ec[i-1] \
                        + (dx/dy+3/4/(dx/dy))*(eps_bot+eps_si)/eps_si*Ec[i] \
                        - 3/8*(dy/dx)*(eps_bot+eps_si)/eps_si*Ec[i+1] \
                        - 1/8*(dy/dx)*eps_bot/eps_si*Ec[i+Nx-1] \
                        - (dx/dy-1/4/(dx/dy))*eps_bot/eps_si*Ec[i+Nx] \
                        - 1/8*(dy/dx)*eps_bot/eps_si*Ec[i+Nx+1]
                elif(i >= (Ntotal-Nx) and i < (Ntotal-Nx+Lsda)):
                    F[i] = Ec[i]-Ec[i-Nx]
                elif(i >= (Ntotal-Nx+Lsda) and i<=(Ntotal-Nx+Lsda+Lg_bota)):
                    F[i] = Ec[i]-Eg2
                elif(i >= (Ntotal-Nx+1+Lsda+Lg_bota) and i < Ntotal):
                    F[i] = Ec[i]-Ec[i-Nx]
                else:
                    F[i] = -(dx/dy)*eps_bot/eps_si*Ec[i-Nx] \
                          - (dy/dx)*eps_bot/eps_si*Ec[i-1] \
                          + 2*(dx/dy+dy/dx)*eps_bot/eps_si*Ec[i] \
                          - (dy/dx)*eps_bot/eps_si*Ec[i+1] \
                          - (dx/dy)*eps_bot/eps_si*Ec[i+Nx]

            #####################Specify the F matrix in the silicon film region################
            for i in range(Nx*(t_topa+1),(Ntotal-Nx*(t_bota+1)+1)-1):
                F[i] = -(dx/dy)*Ec[i-Nx]-(dy/dx)*Ec[i-1]+2*(dx/dy+dy/dx)*Ec[i]+dummy_fun[i]-(dy/dx)*Ec[i+1]-(dx/dy)*Ec[i+Nx]

            #***************`*Modify the F matrix at the right and left boundaries***************
            i_l = 0
            i_r = Nx-1
            for j in range(0,Ny):
                if j == 0:
                    F[i_l]=2*Ec[i_l]-Ec[i_l+1]-Ec[i_l+Nx]
                    F[i_r]=2*Ec[i_r]-Ec[i_r-1]-Ec[i_r+Nx]

                elif(j > 0 and j < (round(Nx*(t_topa+1)/Nx)) - 1):
                    F[i_l]=Ec[i_l]-Ec[i_l+1]
                    F[i_r]=Ec[i_r]-Ec[i_r-1]

                elif(j >= round(Nx*(t_topa)/Nx) and j < round(Ntotal-Nx*t_bota/Nx)):
                    F[i_l]=Ec[i_l]-Ec[i_l+1]
                    F[i_r]=Ec[i_r]-Ec[i_r-1]

                elif(j >= round(Ntotal-Nx*t_bota/Nx) and j < Ny-1):
                    F[i_l]=Ec[i_l]-Ec[i_l+1]
                    F[i_r]=Ec[i_r]-Ec[i_r-1]

                elif(j == Ny-1 and ((Ntotal-Nx+1)<(Ntotal-Nx+1+Lsda)) and ((Ntotal-Nx+1+Lsda+Lg_bota)<Ntotal)):
                    F[i_l]=2*Ec[i_l]-Ec[i_l+1]-Ec[i_l-Nx]
                    F[i_r]=2*Ec[i_r]-Ec[i_r-1]-Ec[i_r-Nx]

                i_l = 1+j*Nx
                i_r = (j+1)*Nx


    ##############################END OF EVALUATING F##################################

    ###############################EVALUATE MF_prime###################################
    # MF_prime matrix in the silicon film region
            for j_row in range (Nx*(t_topa+1)/Nx, (Ntotal-Nx*t_bota)/Nx-2):
                for j_col in range(1,Nx-1):
                    ii = j_row*Nx+j_col
                    MF_prime[ii,ii] = dummy_fun_prime[ii]

    ###############END OF EVALUATION FOR NO PENETRATION INTO THE OXIDE#################

        elif ox_pnt_flag == 1: #(ACCOUNTING FOR ELECTRON PENETRATION INTO OXIDE REGIONS)

    #################################EVALUATE F########################################

    ##########################Top gate insulator region################################
            for i in range(0, Nx*(t_topa+1)):
                if(i >= 0 and i < Lsda+1-1) :
                    F[i] = Ec[i]-Ec[i+Nx]
                elif(i >= Lsda and i < (Lsda+Lg_topa)+1):
                    F[i] = Ec[i]-Eg1
                elif(i >= (Lsda+Lg_topa)+1 and i<Nx):
                    F[i] = Ec[i]-Ec[i+Nx]
                elif(i >= (Nx*t_topa+1) and i <= Nx*(t_topa+1)):
                    F[i] = -1/8*(dy/dx)*eps_top/eps_si*Ec[i-Nx-1] \
                              - (dx/dy-1/4/(dx/dy))*eps_top/eps_si*Ec[i-Nx] \
                              - 1/8*(dy/dx)*eps_top/eps_si*Ec[i-Nx+1] \
                              - 3/8*(dy/dx)*(eps_top+eps_si)/eps_si*Ec[i-1] \
                              + (dx/dy+3/4/(dx/dy))*(eps_top+eps_si)/eps_si*Ec[i] \
                              + dummy_fun[i] \
                              - 3/8*(dy/dx)*(eps_top+eps_si)/eps_si*Ec[i+1] \
                              - 1/8*(dy/dx)*Ec[i+Nx-1] \
                              - (dx/dy-1/4/(dx/dy))*Ec[i+Nx] \
                              - 1/8*(dy/dx)*Ec[i+Nx+1]
                else:
                    F[i] = -(dx/dy)*eps_top/eps_si*Ec[i-Nx] \
                            - (dy/dx)*eps_top/eps_si*Ec[i-1] \
                            + 2*(dx/dy+dy/dx)*eps_top/eps_si*Ec[i]+dummy_fun[i] \
                            - (dy/dx)*eps_top/eps_si*Ec[i+1] \
                            - (dx/dy)*eps_top/eps_si*Ec[i+Nx]

        ############################Bottom gate insulator region###########################
            for i in range((Ntotal-Nx*(t_bota+1)),Ntotal):
                if(i >= (Ntotal-Nx*(t_bota+1)) and i < Ntotal-Nx*t_bota):
                    F[i] = -1/8*(dy/dx)*Ec[i-Nx-1] \
                        - (dx/dy-1/4/(dx/dy))*Ec[i-Nx] \
                        - 1/8*(dy/dx)*Ec[i-Nx+1] \
                        - 3/8*(dy/dx)*(eps_bot+eps_si)/eps_si*Ec[i-1] \
                        + (dx/dy+3/4/(dx/dy))*(eps_bot+eps_si)/eps_si*Ec[i] \
                        + dummy_fun[i] \
                        - 3/8*(dy/dx)*(eps_bot+eps_si)/eps_si*Ec[i+1] \
                        - 1/8*(dy/dx)*eps_bot/eps_si*Ec[i+Nx-1] \
                        - (dx/dy-1/4/(dx/dy))*eps_bot/eps_si*Ec[i+Nx] \
                        - 1/8*(dy/dx)*eps_bot/eps_si*Ec[i+Nx+1]
                elif(i >= (Ntotal-Nx) and i < (Ntotal-Nx+1+Lsda)-1):
                    F[i] = Ec[i]-Ec[i-Nx]
                elif(i >= (Ntotal-Nx+Lsda) and i < (Ntotal-Nx+1+Lsda+Lg_bota)):
                    F[i] = Ec[i]-Eg2
                elif(i >= (Ntotal-Nx+1+Lsda+Lg_bota) and i < Ntotal):
                    F[i] = Ec[i]-Ec[i-Nx]
                else:
                    F[i] = -(dx/dy)*eps_bot/eps_si*Ec[i-Nx] \
                        - (dy/dx)*eps_bot/eps_si*Ec[i-1] \
                        + 2*(dx/dy+dy/dx)*eps_bot/eps_si*Ec[i] + dummy_fun[i] \
                        - (dy/dx)*eps_bot/eps_si*Ec[i+1] \
                        - (dx/dy)*eps_bot/eps_si*Ec[i+Nx]

        ##################Specify the F matrix in the silicon film region#################
            for i in range(Nx*(t_topa+1), (Ntotal-Nx*(t_bota+1)+1)-1):
                F[i] = -(dx/dy)*Ec[i-Nx]-(dy/dx)*Ec[i-1]+2*(dx/dy+dy/dx)*Ec[i]+dummy_fun[i]-(dy/dx)*Ec[i+1]-(dx/dy)*Ec[i+Nx]

        #***************Modify the F matrix at the right and left boundaries**************
            i_l = 0
            i_r = Nx-1
            for j in range(0, Ny):
                if j == 0:
                    F[i_l] = 2*Ec[i_l]-Ec[i_l+1]-Ec[i_l+Nx]
                    F[i_r] = 2*Ec[i_r]-Ec[i_r-1]-Ec[i_r+Nx]

                elif (j > 0 and j < round(Nx*(t_topa+1)/Nx) - 1):
                    F[i_l] = Ec[i_l]-Ec[i_l+1]
                    F[i_r] = Ec[i_r]-Ec[i_r-1]

                elif (j >= round(Nx*(t_topa)/Nx) and j < round(Ntotal-Nx*t_bota/Nx)):
                    F[i_l] = Ec[i_l]-Ec[i_l+1]
                    F[i_r] = Ec[i_r]-Ec[i_r-1]

                elif (j >= round(Ntotal-Nx*t_bota/Nx) and j<Ny-1):
                    F[i_l]=Ec[i_l]-Ec[i_l+1]
                    F[i_r]=Ec[i_r]-Ec[i_r-1]

                elif(j == Ny-1 and ((Ntotal-Nx+1)<(Ntotal-Nx+1+Lsda)) and ((Ntotal-Nx+1+Lsda+Lg_bota)<Ntotal)):
                    F[i_l] = 2*Ec[i_l]-Ec[i_l+1]-Ec[i_l-Nx]
                    F[i_r] = 2*Ec[i_r]-Ec[i_r-1]-Ec[i_r-Nx]

                i_l = 1+j*Nx
                i_r = (j+1)*Nx

        ##########################END OF EVALUATING F###################################

        #############################EVALUATE MF_prime##################################
        # MF_prime matrix in the silicon film region
            for i in range(Nx, (Ntotal-Nx+1)-1):
                MF_prime[i, i] = dummy_fun_prime[i]

            for j in range(1, (Ny-1)):
                MF_prime[(j-1)*Nx+1, (j-1)*Nx+1] = 0
                MF_prime[j*Nx, j*Nx] = 0
        ###############END OF EVALUATION FOR PENETRATION INTO OXIDE ###################
        ##############END OF ELECTRON PENETRATION INTO OXIDE OPTION SWITCH#############

        MF_prime = F_prime+sparse.csr_matrix(MF_prime)
        #######################END OF EVALUATING MF_prime##############################

        ############################SOLVING FOR delta_Ec###############################
        delta_Ec = - np.linalg.solve(sparse.csr_matrix(MF_prime), sparse.csr_matrix(F))

        for i in range(0, Ntotal):
            if abs(delta_Ec(i)) <= 1:
                delta_Ec[i] = delta_Ec[i]
            elif 1<abs(delta_Ec(i)) and abs(delta_Ec(i)) <3.7:
                delta_Ec[i] = np.sign(delta_Ec(i))*np.power(abs(delta_Ec(i)),0.20)
            elif abs(delta_Ec(i)) >= 3.7:
                delta_Ec[i] = np.sign(delta_Ec[i])*np.log(abs(delta_Ec[i]))
        ##########################END OF SOLVING FOR delta_Ec#########################

        Ec = Ec+delta_Ec
        error_inner = max(abs(np.real(F)))
        print '%s %e \n' % ('error_inner = ', error_inner)
        max_delta_Ec = max(abs(np.real(delta_Ec.toarray())))
        MF_prime = np.zeros((Ntotal, Ntotal))
        F = np.zeros((Ntotal, 1))


    ##########################END OF INNER LOOP (WHILE) #######################
    return  Ec
    ###########################################################################
    ###########################END OF FUNCTION POISSON#########################
    ###########################################################################