###########################################################################
#########2D Poisson equation solver using Newton approach (Zhibin Ren 7-18-00)
###########################################################################

from scipy import sparse
from readinput import *
from dummy import dummy


def poisson(spNd, spFn, Ec_old, F_prime, div_avd, charge_fac, Eg1, Eg2, Es, Ed, Nx, Ny, Ntotal):
    Temp = Te
    transport_model = transportmodel.value
    fermi_flag = fermiflag1.value

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
        dummy_fun_prime = charge_fac * dummy_prime((Fn_real-Ec)/(k_B*Temp/q), dummy_flag,fermi_flag)/(k_B*Temp/q)


        ########################THE END OF DUMMY VARIABLE DEFINITION###############

        if ox_pnt_flag == 0: # NO ELECTRON PENETRATION INTO OXIDE REGIONS:

            #################################EVALUATE F################################
            # hpal: modified - F calculated using F_prime
            F = full(F_prime*Ec)

            # Fill in the boundary conditions
            for idex in range (Lsda+1,Lsda+Lg_topa+2):
                F[idex-1] = Ec[idex-1]-Eg1

            for idex in range ((Ntotal-Nx+1+Lsda),(Ntotal-Lsda+1)):
                F[idex-1] = Ec[idex-1]-Eg2

            # Fill in the non-boundary area
            for i in range ((t_topa+1),(t_topa+t_sia)):
                for j in range(2,Nx):
                    idex = Nx*i+j
                    F[idex-1] += dummy_fun(idex)

            ###############################EVALUATE MF_prime###########################
            # MF_prime matrix in the silicon film region
            for j_row in range (Nx*(t_topa+1)/Nx,(Ntotal-Nx*t_bota)/Nx-1):
                for j_col in range (2,Nx):
                    ii = j_row * Nx + j_col
                    MF_prime[ii-1,ii-1] = dummy_fun_prime(ii)

        elif ox_pnt_flag == 1:# ACCOUNTING FOR ELECTRON PENETRATION INTO OXIDE REGIONS

            #################################EVALUATE F################################
            # hpal: modified - F calculated using F_prime
            F = full(F_prime*Ec)

            for idex in range(Lsda+1,Lsda+Lg_topa+2):
                F[idex-1] = Ec[idex-1]-Eg1

            for idex in range ((Ntotal-Nx+1+Lsda),(Ntotal-Lsda+1)):
                F[idex-1] = Ec[idex-1]-Eg2

            for i in range(1,Ny-1):
                for j in range (2,Nx-1):
                    idex = Nx * i + j
                    F[idex-1]=F[idex-1]+dummy_fun[idex-1]

            #############################EVALUATE MF_prime#####################
            #MF_prime matrix in the silicon film region
            for idex in range (Nx+1,(Ntotal-Nx+1)):
                MF_prime[idex-1,idex-1] = dummy_fun_prime[idex-1]

            for jdex in range (2,Ny):
                MF_prime[(jdex-1)*Nx,(jdex-1)*Nx] = 0
                MF_prime[jdex*Nx-1,jdex*Nx-1] = 0

        #######################################################################
        # Debugger ##########
        #######################################################################
        # if debug == 1:
        #     # Plot the contour of F
        #     if debug_counter == 0:
        #         debug_counter += 1
        #         V=np.linspace(1,error_inner/criterion_inner,20)
        #
        #     F_2D = np.reshape(abs(np.real(F)),Nx,Ny)./criterion_inner #Rohit - check reshape compatability
        #     figure(100)
        #     contourf(F_2D, V)
        #     colorbar
        #
        #     pause
        # end

        ##############END OF ELECTRON PENETRATION INTO OXIDE OPTION SWITCH#####

        MF_prime = F_prime+sparse.csr_matrix(MF_prime)

        #######################END OF EVALUATING MF_prime######################

        ############################SOLVING FOR delta_Ec#######################
        # XW - sometimes, MF_prime will become singular due to incorrect Nc
        # value. Adjust Nc to avoid this problem.

        #condition_MF_prime = condest(MF_prime)
        #condition_F = cond(F)


        delta_Ec = -sparse.csr_matrix(MF_prime)/sparse.csr_matrix(F)  #Rohit - \ or / ???

        # nonlinear convergence
        # XW - I do believe this is a limitation on line search to avoid
        # singularity problem
        for idex in range (1,Ntotal+1):
            if abs(delta_Ec[idex-1]) <= 1:
                delta_Ec[idex-1] = delta_Ec[idex-1]
            elif 1 < abs(delta_Ec[idex-1]) and abs(delta_Ec[idex-1]) < 3.7:
                delta_Ec[idex-1] = np.sign(delta_Ec[idex-1])*(abs(delta_Ec[idex-1]).^1/5)  # Rohit modify element wise multiplication from matlabn to python
            elif abs(delta_Ec[idex-1]) >= 3.7:
                delta_Ec[idex-1] = np.sign(delta_Ec[idex-1])*np.log(abs(delta_Ec[idex-1]))
        ##########################END OF SOLVING FOR delta_Ec######################

        Ec = Ec+delta_Ec
        error_inner = max(abs(np.real(F)))
        print '%s %e \n' % ('error_inner = ',error_inner)
        max_delta_Ec = max(abs(np.real(full(delta_Ec))))
        MF_prime = np.zeros(Ntotal,Ntotal)

    return Ec
    ##########################END OF INNER LOOP (WHILE) #######################
    ###########################################################################
    ###########################END OF FUNCTION POISSON#########################
    ###########################################################################